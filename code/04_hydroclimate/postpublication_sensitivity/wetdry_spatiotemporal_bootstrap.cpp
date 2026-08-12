// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix wetdry_spatiotemporal_bootstrap_cpp(
    const IntegerMatrix& wet_count,
    const NumericMatrix& wet_actual,
    const NumericMatrix& wet_iwv,
    const NumericMatrix& wet_z,
    const NumericMatrix& wet_interaction,
    const IntegerMatrix& dry_count,
    const NumericMatrix& dry_actual,
    const NumericMatrix& dry_iwv,
    const NumericMatrix& dry_z,
    const NumericMatrix& dry_interaction,
    const IntegerVector& block_index,
    const NumericVector& area_weight,
    const IntegerMatrix& temporal_year_counts,
    const IntegerMatrix& spatial_block_counts,
    int n_threads = 1) {

  const int n_cell = wet_count.nrow();
  const int n_year = wet_count.ncol();
  const int B = temporal_year_counts.nrow();
  const int n_block = spatial_block_counts.ncol();

  if (dry_count.nrow() != n_cell || dry_count.ncol() != n_year ||
      wet_actual.nrow() != n_cell || wet_actual.ncol() != n_year ||
      wet_iwv.nrow() != n_cell || wet_iwv.ncol() != n_year ||
      wet_z.nrow() != n_cell || wet_z.ncol() != n_year ||
      wet_interaction.nrow() != n_cell ||
      wet_interaction.ncol() != n_year ||
      dry_actual.nrow() != n_cell || dry_actual.ncol() != n_year ||
      dry_iwv.nrow() != n_cell || dry_iwv.ncol() != n_year ||
      dry_z.nrow() != n_cell || dry_z.ncol() != n_year ||
      dry_interaction.nrow() != n_cell ||
      dry_interaction.ncol() != n_year ||
      block_index.size() != n_cell || area_weight.size() != n_cell ||
      temporal_year_counts.ncol() != n_year ||
      spatial_block_counts.nrow() != B) {
    stop("Dimension mismatch in wet/dry bootstrap inputs.");
  }

  NumericMatrix out(B, 8);
  colnames(out) = CharacterVector::create(
    "Wet_actual", "Wet_IWV_channel", "Wet_z_channel",
    "Wet_interaction_channel", "Dry_actual", "Dry_IWV_channel",
    "Dry_z_channel", "Dry_interaction_channel"
  );

#ifdef _OPENMP
  omp_set_num_threads(std::max(1, n_threads));
#pragma omp parallel for schedule(static)
#endif
  for (int b = 0; b < B; ++b) {
    double wet_num[4] = {0.0, 0.0, 0.0, 0.0};
    double dry_num[4] = {0.0, 0.0, 0.0, 0.0};
    double wet_den_space = 0.0;
    double dry_den_space = 0.0;

    for (int i = 0; i < n_cell; ++i) {
      const int block = block_index[i] - 1;
      if (block < 0 || block >= n_block) continue;
      const int spatial_mult = spatial_block_counts(b, block);
      const double aw = area_weight[i];
      if (spatial_mult <= 0 || !R_finite(aw) || aw <= 0.0) continue;
      const double sw = aw * static_cast<double>(spatial_mult);

      int nw = 0;
      int nd = 0;
      double wsum[4] = {0.0, 0.0, 0.0, 0.0};
      double dsum[4] = {0.0, 0.0, 0.0, 0.0};

      for (int y = 0; y < n_year; ++y) {
        const int tm = temporal_year_counts(b, y);
        if (tm <= 0) continue;
        const int cw = wet_count(i, y);
        const int cd = dry_count(i, y);
        nw += tm * cw;
        nd += tm * cd;
        const double mt = static_cast<double>(tm);
        wsum[0] += mt * wet_actual(i, y);
        wsum[1] += mt * wet_iwv(i, y);
        wsum[2] += mt * wet_z(i, y);
        wsum[3] += mt * wet_interaction(i, y);
        dsum[0] += mt * dry_actual(i, y);
        dsum[1] += mt * dry_iwv(i, y);
        dsum[2] += mt * dry_z(i, y);
        dsum[3] += mt * dry_interaction(i, y);
      }

      if (nw > 0) {
        const double inv = 1.0 / static_cast<double>(nw);
        wet_den_space += sw;
        for (int k = 0; k < 4; ++k) wet_num[k] += sw * wsum[k] * inv;
      }
      if (nd > 0) {
        const double inv = 1.0 / static_cast<double>(nd);
        dry_den_space += sw;
        for (int k = 0; k < 4; ++k) dry_num[k] += sw * dsum[k] * inv;
      }
    }

    for (int k = 0; k < 4; ++k) {
      out(b, k) = wet_den_space > 0.0 ? wet_num[k] / wet_den_space : NA_REAL;
      out(b, 4 + k) =
        dry_den_space > 0.0 ? dry_num[k] / dry_den_space : NA_REAL;
    }
  }

  return out;
}
