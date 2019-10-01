// Common definitions for the model block

// -- Priors --
// Study-specific baselines
prior_select_lp(mu, prior_intercept_type, prior_intercept_location, prior_intercept_scale, prior_intercept_df);
// Treatment effects
prior_select_lp(d, prior_trt_type, prior_trt_location, prior_trt_scale, prior_trt_df);
// Regression parameters
prior_select_lp(beta, prior_reg_type, prior_reg_location, prior_reg_scale, prior_reg_df);
// Heterogeneity
prior_select_lp(to_vector(tau), prior_het_type, prior_het_location, prior_het_scale, prior_het_df);

// -- Random effects --
u_delta ~ std_normal();

// -- AgD model (contrast-based) --
if (ni_agd_contrast) {
  agd_contrast_y ~ multi_normal(eta_agd_contrast_bar, agd_contrast_Sigma);
}
