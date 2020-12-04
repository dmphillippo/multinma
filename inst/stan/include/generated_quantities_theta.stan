// Generated quantities for AgD arm integration error
// Used for all univariate outcomes, multivariate outcomes need special attention

// Integration error for AgD arm
for (i in 1:ni_agd_arm) {
  for (j in 1:n_int_thin) {
    theta_bar_cum_agd_arm[(i-1)*n_int_thin + j] = mean(theta_agd_arm_ii[(1 + (i-1)*nint):((i-1)*nint + j*int_thin)]);
  }
}
