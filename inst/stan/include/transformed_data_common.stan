// Common definitions for the transformed data block

// -- Random effects --
// number of random effects
int<lower=0> n_delta = max(which_RE);
// RE MVN mean and correlations
vector[n_delta] RE_mu = rep_vector(0, n_delta);
// Cholesky decomposition of RE MVN correlations
cholesky_factor_corr[n_delta] RE_L = cholesky_decompose(RE_cor);


// Total number of data points
int totni = ni_ipd + nint * (ni_agd_arm + ni_agd_contrast);
int totns = ns_ipd + ns_agd_arm + ns_agd_contrast;

// Split Q matrix into IPD and AgD rows
matrix[ni_ipd, nX] Q_ipd = Q[1:ni_ipd];
matrix[nint * ni_agd_arm, nX] Q_agd_arm = Q[ni_ipd + (1:ni_agd_arm)];
matrix[nint * ni_agd_contrast, nX] Q_agd_contrast = Q[ni_ipd + ni_agd_arm + (1:ni_agd_contrast)];

// Priors on transformed QR scale
// Equivalent to allbeta ~ normal(0, 100)
//vector[totnpar] beta_tilde_mu0 = rep_vector(0, totnpar);
//matrix[totnpar, totnpar] beta_tilde_sigma0 = 100^2 * crossprod(R_inv);

// nint/int_thin for numerical integration checks
int n_int_thin = nint / int_thin;

