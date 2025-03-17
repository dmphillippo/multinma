// Common definitions for the model block

// -- Priors --
// Study-specific baselines
prior_select_lp(mu, prior_intercept_dist, prior_intercept_location, prior_intercept_scale, prior_intercept_df);
// Treatment effects
  if (class_effects == 0) {
    prior_select_lp(d, prior_trt_dist, prior_trt_location, prior_trt_scale, prior_trt_df);
  } else {
    // Priors for class mean parameters
    prior_select_lp(class_mean, prior_class_mean_dist, prior_class_mean_location, prior_class_mean_scale, prior_class_mean_df);

    // Priors for class standard deviation parameters
    prior_select_lp(class_sd, prior_class_sd_dist, prior_class_sd_location, prior_class_sd_scale, prior_class_sd_df);
}

// Regression parameters
prior_select_lp(beta, prior_reg_dist, prior_reg_location, prior_reg_scale, prior_reg_df);

// Node-splitting - just use prior for d
prior_select_lp(omega, prior_trt_dist, prior_trt_location, prior_trt_scale, prior_trt_df);

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

// Class effects model for the d's
if (class_effects) {
  for (t in 1:(nt-1)) {
    if (which_CE[t] == 0) {
      // Priors for treatments without a class effect or in a single-occupancy class
      prior_select2_lp(d[t], prior_trt_dist, prior_trt_location, prior_trt_scale, prior_trt_df);
    } else {
      // Add dummy information to unused allbeta parameters to avoid Stan warnings
      allbeta[(totns + t)] ~ std_normal();
    }
  }
}

// Draw auxiliary variables from standard normal
z_class ~ std_normal();
