// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

static arma::mat raw_columns(const NumericMatrix& I, const NumericMatrix& Z,
                             const NumericMatrix& Int, const NumericMatrix& W,
                             int row, const std::vector<int>& idx, int model) {
  int p = model == 0 ? 1 : (model == 2 ? 3 : 2);
  arma::mat X(idx.size(), p);
  for (size_t a = 0; a < idx.size(); ++a) {
    int t = idx[a];
    if (model == 0) X(a, 0) = I(row, t);
    if (model == 1) { X(a, 0) = I(row, t); X(a, 1) = Z(row, t); }
    if (model == 2) { X(a, 0) = I(row, t); X(a, 1) = Z(row, t); X(a, 2) = Int(row, t); }
    if (model == 3) { X(a, 0) = I(row, t); X(a, 1) = W(row, t); }
  }
  return X;
}

static void standardize_training(arma::mat& tr, arma::mat& te) {
  for (arma::uword j = 0; j < tr.n_cols; ++j) {
    double mu = arma::mean(tr.col(j));
    double sd = arma::stddev(tr.col(j));
    if (!(sd > 0) || !std::isfinite(sd)) sd = 1.0;
    tr.col(j) = (tr.col(j) - mu) / sd;
    te.col(j) = (te.col(j) - mu) / sd;
  }
}

// [[Rcpp::export]]
List cv_oof_predictions_training_only_cpp(NumericMatrix Y, NumericMatrix I,
                                          NumericMatrix Z, NumericMatrix Int,
                                          NumericMatrix W, IntegerVector fold,
                                          int gap = 3) {
  const int nr = Y.nrow(), nt = Y.ncol(), nf = max(fold);
  NumericMatrix p0(nr, nt), p1(nr, nt), p2(nr, nt), p3(nr, nt);
  std::fill(p0.begin(), p0.end(), NA_REAL);
  std::fill(p1.begin(), p1.end(), NA_REAL);
  std::fill(p2.begin(), p2.end(), NA_REAL);
  std::fill(p3.begin(), p3.end(), NA_REAL);

  #pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < nr; ++i) {
    std::vector<int> valid;
    valid.reserve(nt);
    for (int t = 0; t < nt; ++t) {
      if (R_finite(Y(i,t)) && R_finite(I(i,t)) && R_finite(Z(i,t)) &&
          R_finite(Int(i,t)) && R_finite(W(i,t))) valid.push_back(t);
    }
    if (valid.size() < 36) continue;

    for (int ff = 1; ff <= nf; ++ff) {
      int lo = nt, hi = -1;
      for (int t = 0; t < nt; ++t) if (fold[t] == ff) {
        lo = std::min(lo, t); hi = std::max(hi, t);
      }
      if (hi < 0) continue;
      std::vector<int> tr, te;
      for (int t : valid) {
        if (fold[t] == ff) te.push_back(t);
        else if (t < lo-gap || t > hi+gap) tr.push_back(t);
      }
      if (te.empty() || tr.size() < 20) continue;

      arma::vec ytr(tr.size());
      for (size_t a = 0; a < tr.size(); ++a) ytr[a] = Y(i, tr[a]);
      double ymu = arma::mean(ytr), ysd = arma::stddev(ytr);
      if (!(ysd > 0) || !std::isfinite(ysd)) continue;
      arma::vec ystd = (ytr - ymu) / ysd;

      for (int m = 0; m < 4; ++m) {
        arma::mat Xtr = raw_columns(I, Z, Int, W, i, tr, m);
        arma::mat Xte = raw_columns(I, Z, Int, W, i, te, m);
        standardize_training(Xtr, Xte);
        if (m == 3) {
          double bwi = arma::dot(Xtr.col(0), Xtr.col(1)) /
                       std::max(1e-15, arma::dot(Xtr.col(0), Xtr.col(0)));
          Xtr.col(1) -= bwi * Xtr.col(0);
          Xte.col(1) -= bwi * Xte.col(0);
        }
        Xtr.insert_cols(0, arma::vec(Xtr.n_rows, arma::fill::ones));
        Xte.insert_cols(0, arma::vec(Xte.n_rows, arma::fill::ones));
        arma::vec beta;
        bool ok = arma::solve(beta, Xtr, ystd);
        if (!ok || beta.n_elem != Xtr.n_cols) beta = arma::pinv(Xtr) * ystd;
        if (beta.n_elem != Xtr.n_cols) continue;
        arma::vec phat = (Xte * beta) * ysd + ymu;
        for (size_t a = 0; a < te.size(); ++a) {
          int t = te[a];
          if (m == 0) p0(i,t) = phat[a];
          if (m == 1) p1(i,t) = phat[a];
          if (m == 2) p2(i,t) = phat[a];
          if (m == 3) p3(i,t) = phat[a];
        }
      }
    }
  }
  return List::create(_["observed"] = Y, _["pred_M0"] = p0,
                      _["pred_M1"] = p1, _["pred_M2"] = p2,
                      _["pred_M3"] = p3);
}
