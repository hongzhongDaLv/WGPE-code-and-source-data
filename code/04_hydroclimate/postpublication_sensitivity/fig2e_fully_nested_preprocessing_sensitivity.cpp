// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <limits>

using namespace Rcpp;

struct Processed {
  arma::vec tr;
  arma::vec te;
  bool ok;
  int fail_code; // 1 climatology, 2 trend, 3 non-finite
};

static Processed preprocess_fold(
    const NumericMatrix& x, int row,
    const std::vector<int>& tr_idx,
    const std::vector<int>& te_idx,
    const IntegerVector& month,
    const NumericVector& tdec,
    bool remove_calendar_climatology) {

  Processed out;
  out.ok = false;
  out.fail_code = 0;
  out.tr.set_size(tr_idx.size());
  out.te.set_size(te_idx.size());

  double clim[12];
  int nclim[12];
  for (int m = 0; m < 12; ++m) {
    clim[m] = 0.0;
    nclim[m] = 0;
  }

  if (remove_calendar_climatology) {
    for (int t : tr_idx) {
      const int m = month[t] - 1;
      if (m < 0 || m >= 12 || !R_finite(x(row, t))) {
        out.fail_code = 1;
        return out;
      }
      clim[m] += x(row, t);
      nclim[m]++;
    }
    for (int m = 0; m < 12; ++m) {
      if (nclim[m] < 1) {
        out.fail_code = 1;
        return out;
      }
      clim[m] /= static_cast<double>(nclim[m]);
    }
  } else {
    for (int m = 0; m < 12; ++m) clim[m] = 0.0;
  }

  double mt = 0.0, my = 0.0;
  for (int t : tr_idx) {
    const double yy = x(row, t) - clim[month[t] - 1];
    mt += tdec[t];
    my += yy;
  }
  mt /= static_cast<double>(tr_idx.size());
  my /= static_cast<double>(tr_idx.size());

  double sxx = 0.0, sxy = 0.0;
  for (int t : tr_idx) {
    const double dt = tdec[t] - mt;
    const double yy = x(row, t) - clim[month[t] - 1];
    sxx += dt * dt;
    sxy += dt * (yy - my);
  }
  if (!(sxx > 0.0) || !std::isfinite(sxx) || !std::isfinite(sxy)) {
    out.fail_code = 2;
    return out;
  }
  const double slope = sxy / sxx;
  const double intercept = my - slope * mt;

  for (size_t a = 0; a < tr_idx.size(); ++a) {
    const int t = tr_idx[a];
    out.tr[a] = x(row, t) - clim[month[t] - 1] -
      (intercept + slope * tdec[t]);
  }
  for (size_t a = 0; a < te_idx.size(); ++a) {
    const int t = te_idx[a];
    out.te[a] = x(row, t) - clim[month[t] - 1] -
      (intercept + slope * tdec[t]);
  }
  if (!out.tr.is_finite() || !out.te.is_finite()) {
    out.fail_code = 3;
    return out;
  }
  out.ok = true;
  return out;
}

static bool standardize_training(arma::vec& tr, arma::vec& te,
                                 double& mu, double& sd) {
  mu = arma::mean(tr);
  arma::vec centered = tr - mu;
  const double denom = std::max(1.0, static_cast<double>(tr.n_elem - 1));
  sd = std::sqrt(arma::dot(centered, centered) / denom);
  if (!(sd > 1e-12) || !std::isfinite(sd)) return false;
  tr = centered / sd;
  te = (te - mu) / sd;
  return tr.is_finite() && te.is_finite();
}

static bool fit_predict(const arma::mat& xtr_raw,
                        const arma::mat& xte_raw,
                        const arma::vec& ytr_std,
                        double ymu, double ysd,
                        arma::vec& pred_physical) {
  arma::mat xtr(xtr_raw.n_rows, xtr_raw.n_cols + 1);
  arma::mat xte(xte_raw.n_rows, xte_raw.n_cols + 1);
  xtr.col(0).ones();
  xte.col(0).ones();
  xtr.cols(1, xtr.n_cols - 1) = xtr_raw;
  xte.cols(1, xte.n_cols - 1) = xte_raw;
  arma::vec beta;
  bool solved = false;
  try {
    solved = arma::solve(beta, xtr, ytr_std, arma::solve_opts::fast);
  } catch (...) {
    solved = false;
  }
  if (!solved || beta.n_elem != xtr.n_cols || !beta.is_finite()) {
    try {
      beta = arma::pinv(xtr) * ytr_std;
    } catch (...) {
      return false;
    }
  }
  if (beta.n_elem != xtr.n_cols || !beta.is_finite()) return false;
  pred_physical = (xte * beta) * ysd + ymu;
  return pred_physical.is_finite();
}

// [[Rcpp::export]]
List fig2e_fully_nested_oof_cpp(
    NumericMatrix Y,
    NumericMatrix I,
    NumericMatrix Z,
    NumericMatrix W,
    IntegerVector month,
    NumericVector tdec,
    IntegerVector fold,
    bool response_remove_calendar_climatology,
    int gap = 3,
    int min_train = 36) {

  const int nr = Y.nrow();
  const int nt = Y.ncol();
  const int nf = max(fold);
  if (I.nrow() != nr || Z.nrow() != nr || W.nrow() != nr ||
      I.ncol() != nt || Z.ncol() != nt || W.ncol() != nt ||
      month.size() != nt || tdec.size() != nt || fold.size() != nt) {
    stop("Dimension mismatch in fully nested OOF inputs.");
  }

  NumericMatrix r2(nr, 4), rmse(nr, 4);
  NumericMatrix fold_delta(nr, nf * 3), fold_n_test(nr, nf);
  IntegerVector n_valid_raw(nr), n_oof_common(nr), n_oof_M0(nr),
    n_oof_M1(nr), n_oof_M2(nr), n_oof_M3(nr), duplicate_oof(nr),
    folds_attempted(nr), folds_success(nr), fail_insufficient(nr),
    fail_preprocessing(nr), fail_zero_variance(nr), fail_model(nr),
    min_train_success(nr), min_test_success(nr);

  std::fill(r2.begin(), r2.end(), NA_REAL);
  std::fill(rmse.begin(), rmse.end(), NA_REAL);
  std::fill(fold_delta.begin(), fold_delta.end(), NA_REAL);
  std::fill(fold_n_test.begin(), fold_n_test.end(), NA_REAL);

  #pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < nr; ++i) {
    std::vector<int> valid;
    valid.reserve(nt);
    for (int t = 0; t < nt; ++t) {
      if (R_finite(Y(i, t)) && R_finite(I(i, t)) &&
          R_finite(Z(i, t)) && R_finite(W(i, t))) {
        valid.push_back(t);
      }
    }
    n_valid_raw[i] = valid.size();
    std::vector<int> seen(nt, 0);
    double sum_y = 0.0, sum_y2 = 0.0;
    double sse[4] = {0.0, 0.0, 0.0, 0.0};
    int nobs = 0;
    int min_tr = nt + 1, min_te = nt + 1;

    for (int ff = 1; ff <= nf; ++ff) {
      int lo = nt, hi = -1;
      for (int t = 0; t < nt; ++t) {
        if (fold[t] == ff) {
          lo = std::min(lo, t);
          hi = std::max(hi, t);
        }
      }
      if (hi < 0) continue;

      std::vector<int> tr, te;
      tr.reserve(valid.size());
      te.reserve(64);
      for (int t : valid) {
        if (fold[t] == ff) {
          te.push_back(t);
        } else if (t < lo - gap || t > hi + gap) {
          tr.push_back(t);
        }
      }
      folds_attempted[i]++;
      if (te.empty() || static_cast<int>(tr.size()) < min_train) {
        fail_insufficient[i]++;
        continue;
      }

      Processed py = preprocess_fold(Y, i, tr, te, month, tdec,
                                     response_remove_calendar_climatology);
      Processed pi = preprocess_fold(I, i, tr, te, month, tdec, true);
      Processed pz = preprocess_fold(Z, i, tr, te, month, tdec, true);
      Processed pw = preprocess_fold(W, i, tr, te, month, tdec, true);
      if (!py.ok || !pi.ok || !pz.ok || !pw.ok) {
        fail_preprocessing[i]++;
        continue;
      }

      double ymu, ysd, imu, isd, zmu, zsd, wmu, wsd;
      arma::vec ytr_std = py.tr, yte_std = py.te;
      arma::vec itr = pi.tr, ite = pi.te;
      arma::vec ztr = pz.tr, zte = pz.te;
      arma::vec wtr = pw.tr, wte = pw.te;
      if (!standardize_training(ytr_std, yte_std, ymu, ysd) ||
          !standardize_training(itr, ite, imu, isd) ||
          !standardize_training(ztr, zte, zmu, zsd) ||
          !standardize_training(wtr, wte, wmu, wsd)) {
        fail_zero_variance[i]++;
        continue;
      }

      arma::vec int_tr = itr % ztr;
      arma::vec int_te = ite % zte;
      double int_mu, int_sd;
      if (!standardize_training(int_tr, int_te, int_mu, int_sd)) {
        fail_zero_variance[i]++;
        continue;
      }

      const double denom_i = arma::dot(itr, itr);
      if (!(denom_i > 1e-12) || !std::isfinite(denom_i)) {
        fail_zero_variance[i]++;
        continue;
      }
      const double b_wi = arma::dot(itr, wtr) / denom_i;
      arma::vec wres_tr = wtr - b_wi * itr;
      arma::vec wres_te = wte - b_wi * ite;
      if (!wres_tr.is_finite() || !wres_te.is_finite()) {
        fail_zero_variance[i]++;
        continue;
      }

      arma::mat xtr[4], xte[4];
      xtr[0].set_size(tr.size(), 1);
      xte[0].set_size(te.size(), 1);
      xtr[0].col(0) = itr; xte[0].col(0) = ite;
      xtr[1].set_size(tr.size(), 2);
      xte[1].set_size(te.size(), 2);
      xtr[1].col(0) = itr; xtr[1].col(1) = ztr;
      xte[1].col(0) = ite; xte[1].col(1) = zte;
      xtr[2].set_size(tr.size(), 3);
      xte[2].set_size(te.size(), 3);
      xtr[2].col(0) = itr; xtr[2].col(1) = ztr; xtr[2].col(2) = int_tr;
      xte[2].col(0) = ite; xte[2].col(1) = zte; xte[2].col(2) = int_te;
      xtr[3].set_size(tr.size(), 2);
      xte[3].set_size(te.size(), 2);
      xtr[3].col(0) = itr; xtr[3].col(1) = wres_tr;
      xte[3].col(0) = ite; xte[3].col(1) = wres_te;

      arma::vec pred[4];
      bool all_models_ok = true;
      for (int m = 0; m < 4; ++m) {
        if (!fit_predict(xtr[m], xte[m], ytr_std, ymu, ysd, pred[m])) {
          all_models_ok = false;
          break;
        }
      }
      if (!all_models_ok) {
        fail_model[i]++;
        continue;
      }

      const double fold_ymean = arma::mean(py.te);
      double fold_sst = 0.0, fold_sse[4] = {0.0, 0.0, 0.0, 0.0};
      for (size_t a = 0; a < te.size(); ++a) {
        const int t = te[a];
        seen[t]++;
        const double yy = py.te[a];
        sum_y += yy;
        sum_y2 += yy * yy;
        nobs++;
        const double dc = yy - fold_ymean;
        fold_sst += dc * dc;
        for (int m = 0; m < 4; ++m) {
          const double err = yy - pred[m][a];
          sse[m] += err * err;
          fold_sse[m] += err * err;
        }
      }
      if (fold_sst > 1e-20 && std::isfinite(fold_sst)) {
        fold_delta(i, (ff - 1) * 3 + 0) =
          (fold_sse[0] - fold_sse[1]) / fold_sst;
        fold_delta(i, (ff - 1) * 3 + 1) =
          (fold_sse[0] - fold_sse[2]) / fold_sst;
        fold_delta(i, (ff - 1) * 3 + 2) =
          (fold_sse[0] - fold_sse[3]) / fold_sst;
      }
      fold_n_test(i, ff - 1) = te.size();
      folds_success[i]++;
      min_tr = std::min(min_tr, static_cast<int>(tr.size()));
      min_te = std::min(min_te, static_cast<int>(te.size()));
    }

    int dup = 0, predicted = 0;
    for (int t = 0; t < nt; ++t) {
      if (seen[t] > 1) dup += seen[t] - 1;
      if (seen[t] == 1) predicted++;
    }
    duplicate_oof[i] = dup;
    n_oof_common[i] = predicted;
    n_oof_M0[i] = predicted;
    n_oof_M1[i] = predicted;
    n_oof_M2[i] = predicted;
    n_oof_M3[i] = predicted;
    min_train_success[i] = (min_tr == nt + 1) ? NA_INTEGER : min_tr;
    min_test_success[i] = (min_te == nt + 1) ? NA_INTEGER : min_te;

    if (nobs >= 2) {
      const double sst = sum_y2 - sum_y * sum_y / static_cast<double>(nobs);
      if (sst > 1e-20 && std::isfinite(sst)) {
        for (int m = 0; m < 4; ++m) {
          r2(i, m) = 1.0 - sse[m] / sst;
          rmse(i, m) = std::sqrt(sse[m] / static_cast<double>(nobs));
        }
      }
    }
  }

  return List::create(
    _["r2"] = r2,
    _["rmse"] = rmse,
    _["fold_delta"] = fold_delta,
    _["fold_n_test"] = fold_n_test,
    _["n_valid_raw"] = n_valid_raw,
    _["n_oof_common"] = n_oof_common,
    _["n_oof_M0"] = n_oof_M0,
    _["n_oof_M1"] = n_oof_M1,
    _["n_oof_M2"] = n_oof_M2,
    _["n_oof_M3"] = n_oof_M3,
    _["duplicate_oof"] = duplicate_oof,
    _["folds_attempted"] = folds_attempted,
    _["folds_success"] = folds_success,
    _["fail_insufficient"] = fail_insufficient,
    _["fail_preprocessing"] = fail_preprocessing,
    _["fail_zero_variance"] = fail_zero_variance,
    _["fail_model"] = fail_model,
    _["min_train_success"] = min_train_success,
    _["min_test_success"] = min_test_success
  );
}
