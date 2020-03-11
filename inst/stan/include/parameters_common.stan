// Common definitions for the parameters block
// Parameters (on QR scale if QR = 1)
vector[nX] beta_tilde;

// (Uncorrelated) random effects and heterogeneity SD
vector[n_delta] u_delta;
vector<lower = 0>[RE ? 1 : 0] tau;
