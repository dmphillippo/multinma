// Common definitions for the generated quantities block

// -- Log likelihood and residual deviance calculation --
vector[ni_ipd + ni_agd_arm + ni_agd_contrast] log_lik;
vector[ni_ipd + ni_agd_arm + ni_agd_contrast] resdev;

// -- Estimate integration error --
vector[(ni_agd_arm + ni_agd_contrast) * n_int_thin] theta_bar_cum;

// -- RE shrunken estimate delta --
// Note: These are the individual-level trial-specific treatment effects
vector[RE ? narm_ipd + ni_agd_arm + ni_agd_contrast : 0] delta;

if (RE) {
  for (i in 1:(narm_ipd + ni_agd_arm + ni_agd_contrast)) {
    delta[i] = (trt[i] > 1 ? gamma[trt[i] - 1] : 0) + (which_RE[i] ? f_delta[which_RE[i]] : 0);
  }
  for (i in 1:ni_agd_contrast) {
    if (which_RE[narm_ipd + ni_agd_arm + i] && agd_contrast_trt_b[i] > 1)
      delta[narm_ipd + ni_agd_arm + i] -= gamma[agd_contrast_trt_b[i] - 1];
  }
}

// Integration error for AgD
for (i in 1:ni_agd_arm) {
  for (j in 1:n_int_thin) {
    theta_bar_cum[(i-1)*n_int_thin + j] = mean(theta_agd_arm_ii[(1 + (i-1)*nint):((i-1)*nint + j*int_thin)]);
  }
}
for (i in 1:ni_agd_contrast) {
  for (j in 1:n_int_thin) {
    theta_bar_cum[ni_agd_arm + (i-1)*n_int_thin + j] = mean(eta_agd_contrast_ii[(1 + (i-1)*nint):((i-1)*nint + j*int_thin)]);
  }
}
