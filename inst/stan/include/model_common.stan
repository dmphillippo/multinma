// Common definitions for the model block

// -- Priors --
// Study-specific baselines
mu ~ prior_select(prior_intercept_type, prior_intercept_location, prior_intercept_scale, prior_intercept_df);
// Treatment effects
gamma ~ prior_select(prior_trt_type, prior_trt_location, prior_trt_scale, prior_trt_df);
// Regression parameters
beta ~ prior_select(prior_reg_type, prior_reg_location, prior_reg_scale, prior_reg_df);
// Heterogeneity
tau ~ prior_select(prior_het_type, prior_het_location, prior_het_scale, prior_het_df);

// -- Random effects --
u_delta ~ std_normal();

// -- AgD model (contrast-based) --
agd_contrast_y ~ normal(eta_agd_contrast, agd_contrast_se);
