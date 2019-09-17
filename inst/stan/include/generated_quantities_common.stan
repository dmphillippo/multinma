// Common definitions for the generated quantities block

// -- Log likelihood and residual deviance calculation --
vector[ni_ipd + ni_agd_arm + ni_agd_contrast] log_lik;
vector[ni_ipd + ni_agd_arm + ni_agd_contrast] resdev;

// -- Estimate integration error --
vector[(ni_agd_arm + ni_agd_contrast) * n_int_thin] theta_bar_cum;

// -- RE shrunken estimate delta --
// Note: These are the individual-level trial-specific treatment effects
vector[n_delta ? narm_ipd + ni_agd_arm + ni_agd_contrast : 0] delta;

for (i in 1:(narm_ipd + ni_agd_arm + ni_agd_contrast)) {
  delta[i] = which_re[i] ? 0 : gamma[trt[i] - 1] + f_delta[which_RE[i]];
}

// Integration error for AgD
for (i in 1:ni_agd_arm) {
  for (j in 1:n_int_thin) {
    theta_bar_cum[(i-1)*n_int_thin + j] = mean(theta_ii[(1 + (i-1)*nint):((i-1)*nint + j*int_thin)]);
  }
}
for (i in (ni_agd_arm + 1):(ni_agd_arm + ni_agd_contrast)) {
  for (j in 1:n_int_thin) {
    theta_bar_cum[(i-1)*n_int_thin + j] = mean(eta_ii[(1 + (i-1)*nint):((i-1)*nint + j*int_thin)]);
  }
}
