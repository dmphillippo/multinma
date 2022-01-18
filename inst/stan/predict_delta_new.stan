data {
  int<lower=1> nt; // number of treatments
  corr_matrix[nt - 1] RE_cor; // RE correlation matrix
}
transformed data {
  cholesky_factor_corr[nt - 1] RE_L = cholesky_decompose(RE_cor);
}
parameters {
  vector[nt - 1] d; // treatment effects
  real<lower = 0> tau; // RE heterogeneity sd
}
generated quantities {
  // predicted treatment effect in new study
  vector[nt - 1] delta_new = multi_normal_cholesky_rng(d, tau * RE_L);
}
