// Common definitions for the parameters block
// Parameters (on QR scale if QR = 1)
vector[nX] beta_tilde;

// (Uncorrelated) random effects and heterogeneity SD
vector[n_delta] u_delta;
vector<lower = 0>[RE ? 1 : 0] tau;

// New vector parameters for class effects model
vector[class_effects ? max(which_ce_num) : 0] class_mean; // Class means (zero-length if no exchangeable class effects)
vector<lower=0>[class_effects ? max(which_ce_sd_num) : 0] class_sd; // Class standard deviations (zero-length if no exchangeable class effects)
