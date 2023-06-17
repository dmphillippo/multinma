// Transformed predictors for IPD and AgD arm
// Used for all univariate outcomes, multivariate outcomes need special attention
vector[ni_ipd] theta_ipd; // IPD transformed predictor
vector[nint_max > 1 ? nint * ni_agd_arm : 0] theta_agd_arm_ii; // AgD arm transformed predictor integration points
vector[ni_agd_arm] theta_agd_arm_bar; // AgD arm transformed predictor
