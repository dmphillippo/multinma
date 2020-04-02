// Common definitions for the model block

// -- Priors --
// Study-specific baselines
prior_select_lp(mu, prior_intercept_dist, prior_intercept_location, prior_intercept_scale, prior_intercept_df);
// Treatment effects
prior_select_lp(d, prior_trt_dist, prior_trt_location, prior_trt_scale, prior_trt_df);
// Regression parameters
prior_select_lp(beta, prior_reg_dist, prior_reg_location, prior_reg_scale, prior_reg_df);

// Heterogeneity
if (RE) {
  if (prior_het_type == 1) { // Prior on sd = tau
    prior_select_lp(tau, prior_het_dist, prior_het_location, prior_het_scale, prior_het_df);

  } else {{ // Prior on transformed parameter, with Jacobian adjustment
    // Transform of tau (e.g. to var or prec)
    vector[RE ? 1 : 0] tau_t;

    if (prior_het_type == 2) { // Prior on var = tau^2
      tau_t[1] = tau[1]^2;
      target += log(tau);
    }
    else if (prior_het_type == 3) { // Prior on prec = tau^-2
      tau_t[1] = tau[1]^-2;
      target += -3*log(tau);
    }
    prior_select_lp(tau_t, prior_het_dist, prior_het_location, prior_het_scale, prior_het_df);
  }}
}

// -- Random effects --
u_delta ~ std_normal();

// -- AgD model (contrast-based) --
if (ni_agd_contrast) {
  agd_contrast_y ~ multi_normal(eta_agd_contrast_bar, agd_contrast_Sigma);
}
