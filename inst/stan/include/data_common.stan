// Common definitions for the data block
// -- Constants --
int<lower=0> ns_ipd; // number of IPD studies
int<lower=0> ns_agd_arm; // number of AgD (arm-based) studies
int<lower=0> ns_agd_contrast; // number of AgD (contrast-based) studies
int<lower=0> ns_agd_regression; // number of AgD (regression coefficient) studies

int<lower=0> ni_ipd; // total number of IPD individuals
int<lower=0> ni_agd_arm; // total number of AgD (arm-based) data points
int<lower=0> ni_agd_contrast; // total number of AgD (contrast-based) data points
int<lower=0> ni_agd_regression; // total number of AgD (regression coefficient) data points

// Treatment IDs
int<lower=0> narm_ipd; // Number of IPD arms
array[ni_ipd] int<lower=1> ipd_arm; // Arm indicator for IPD (i.e. picking element of which_RE)
array[narm_ipd] int<lower=1> ipd_trt;
int<lower=0> narm_agd_arm;
array[narm_agd_arm] int<lower=1> agd_arm_trt;
array[ni_agd_contrast] int<lower=1> agd_contrast_trt;
array[ni_agd_contrast] int<lower=1> agd_contrast_trt_b;

// Study IDs
// array[max(ipd_arm)] int<lower=1> ipd_study;
// array[ni_agd_arm] int<lower=1> agd_arm_study;
// array[ni_agd_contrast] int<lower=1> agd_contrast_study;

int<lower=1> nt; // number of treatments
int<lower=0> nX; // number of columns of design matrix

// Integration
int<lower=1> nchains;
int<lower=1,upper=nchains> CHAIN_ID;
array[nchains] int<lower=1> nint_vec; // number of samples for numerical integration (1 = no integration)
int<lower=1> nint_max;
int<lower=0> int_thin; // thinning factor for saved integration points

int<lower=1> link; // link function

// -- AgD (contrast-based) outcomes
vector[ni_agd_contrast] agd_contrast_y;
cov_matrix[ni_agd_contrast ? ni_agd_contrast : 1] agd_contrast_Sigma;

// -- AgD regression coefficients --
vector[ni_agd_regression] agd_regression_est;
cov_matrix[ni_agd_regression ? ni_agd_regression : 1] agd_regression_cov;

// -- Design matrix or thin QR decomposition --
int<lower=0, upper=1> QR; // use QR decomposition (yes = 1)
matrix[ni_ipd + nint_max * (ni_agd_arm + ni_agd_contrast) + ni_agd_regression, nX] X; // X is Q from QR decomposition if QR = 1
matrix[QR ? nX : 0, QR ? nX : 0] R_inv;

// -- Offsets --
int<lower=0, upper=1> has_offset; // Offset flag (yes = 1)
vector[has_offset ? ni_ipd + nint_max * (ni_agd_arm + ni_agd_contrast) : 0] offsets; // Vector of offsets

// -- Random effects --
int<lower=0, upper=1> RE; // Random effects flag (yes = 1)
array[RE ? narm_ipd + narm_agd_arm + ni_agd_contrast + ni_agd_regression : 0] int<lower=0> which_RE; // ID of RE delta for each arm (0 for no RE delta)
corr_matrix[RE ? max(which_RE) : 1] RE_cor; // RE correlation matrix

// -- Node-splitting --
int<lower=0, upper=1> nodesplit; // Node-splitting flag (yes = 1)

// -- Priors --
int<lower=0,upper=3> prior_intercept_dist;
real prior_intercept_location;
real<lower=0> prior_intercept_scale;
real<lower=0> prior_intercept_df;

int<lower=0,upper=3> prior_trt_dist;
real prior_trt_location;
real<lower=0> prior_trt_scale;
real<lower=0> prior_trt_df;

int<lower=0,upper=6> prior_het_dist;
int<lower=1,upper=3> prior_het_type;
real prior_het_location;
real<lower=0> prior_het_scale;
real<lower=0> prior_het_df;

int<lower=0,upper=3> prior_reg_dist;
real prior_reg_location;
real<lower=0> prior_reg_scale;
real<lower=0> prior_reg_df;
