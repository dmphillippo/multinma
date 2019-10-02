// Common definitions for the transformed data block

// -- Random effects --
// number of random effects
int<lower=0> n_delta = RE ? max(which_RE) : 0;
// RE MVN mean and correlations
vector[n_delta] RE_mu = rep_vector(0, n_delta);
// Cholesky decomposition of RE MVN correlations
matrix[0, 0] REdummy; // Use dummy zero-dim matrix to set RE_L to [0x0] when n_delta = 0
cholesky_factor_corr[n_delta] RE_L = n_delta ? cholesky_decompose(RE_cor) : REdummy;

// Total number of data points
int totni = ni_ipd + nint * (ni_agd_arm + ni_agd_contrast);

// Total number of study intercepts (none for contrast-based data)
int totns = ns_ipd + ns_agd_arm; // + ns_agd_contrast;

// Number of IPD arms
// int<lower=0> narm_ipd = ni_ipd ? max(ipd_arm) : 0;

// All treatments vector
int<lower=1> trt[narm_ipd + ni_agd_arm + ni_agd_contrast] = append_array(append_array(ipd_trt, agd_arm_trt), agd_contrast_trt);

// Split Q matrix or X matrix into IPD and AgD rows
matrix[0, nX] Xdummy;
matrix[ni_ipd, nX] X_ipd = ni_ipd ? X[1:ni_ipd] : Xdummy;
matrix[nint * ni_agd_arm, nX] X_agd_arm = ni_agd_arm ? X[(ni_ipd + 1):(ni_ipd + nint * ni_agd_arm)] : Xdummy;
matrix[nint * ni_agd_contrast, nX] X_agd_contrast = ni_agd_contrast ? X[(ni_ipd + nint * ni_agd_arm + 1):(ni_ipd + nint * (ni_agd_arm + ni_agd_contrast))] : Xdummy;

// nint/int_thin for numerical integration checks
int n_int_thin = nint / int_thin;

// Inverse covariance matrix for contrasts
matrix[ns_agd_contrast, ns_agd_contrast] inv_Sigma = inverse_spd(agd_contrast_Sigma);

// Construct number of contrasts in each study for contrast-based AgD by looking at Sigma covariance matrix
// NOTE: Sigma must be block diagonal (i.e. all contrasts for a single study together)
int nc_agd_contrast[ns_agd_contrast];
if (ns_agd_contrast) {
  int s = 1;
  int c = 1;
  for (i in 1:(ni_agd_contrast - 1)) {
    if (agd_contrast_Sigma[i, i+1] == 0) {
      nc_agd_contrast[s] = c;
      s += 1;
      c = 1;
    } else {
      c += 1;
    }
  }
  // for i = ni_agd_contrast
  nc_agd_contrast[s] = c;
}
