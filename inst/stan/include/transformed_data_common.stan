// Common definitions for the transformed data block

// -- Integration --
// Set nint by chain
int<lower=1,upper=nint_max> nint = nint_vec[CHAIN_ID];

// -- Random effects --
// number of random effects
int<lower=0> n_delta = RE ? max(which_RE) : 0;
// RE MVN mean and correlations
vector[n_delta] RE_mu = rep_vector(0, n_delta);
// Cholesky decomposition of RE MVN correlations
matrix[0, 0] REdummy; // Use dummy zero-dim matrix to set RE_L to [0x0] when n_delta = 0
cholesky_factor_corr[n_delta] RE_L = n_delta ? cholesky_decompose(RE_cor) : REdummy;
// Sparse representation
vector[0] wdummy;
array[0] int vudummy;
int RE_L_nz = count_nonzero(RE_L); // Number of non-zero entries
int RE_sparse = RE_L_nz * 1.0 / num_elements(RE_L) <= 0.1; // Use sparse representation? (yes = 1)
vector[RE_sparse ? RE_L_nz : 0] RE_L_w = RE_sparse ? csr_extract_w(RE_L): wdummy; // Non-zero entries
array[RE_sparse ? RE_L_nz : 0] int RE_L_v = RE_sparse ? csr_extract_v(RE_L): vudummy; // v sparse component
array[RE_sparse ? n_delta + 1 : 0] int RE_L_u = RE_sparse ? csr_extract_u(RE_L) : vudummy; // u sparse component

// Total number of data points
// int totni = ni_ipd + nint * (ni_agd_arm + ni_agd_contrast);

// Total number of study intercepts (none for contrast-based data)
int totns = ns_ipd + ns_agd_arm; // + ns_agd_contrast;

// Number of IPD arms
// int<lower=0> narm_ipd = ni_ipd ? max(ipd_arm) : 0;

// All treatments vector
array[narm_ipd + narm_agd_arm + ni_agd_contrast] int<lower=1> trt = append_array(append_array(ipd_trt, agd_arm_trt), agd_contrast_trt);

// Split Q matrix or X matrix into IPD and AgD rows
matrix[0, nX] Xdummy;
matrix[ni_ipd, nX] X_ipd = ni_ipd ? X[1:ni_ipd] : Xdummy;
matrix[nint_max * ni_agd_arm, nX] X_agd_arm = ni_agd_arm ? X[(ni_ipd + 1):(ni_ipd + nint_max * ni_agd_arm)] : Xdummy;
matrix[nint_max * ni_agd_contrast, nX] X_agd_contrast = ni_agd_contrast ? X[(ni_ipd + nint_max * ni_agd_arm + 1):(ni_ipd + nint_max * (ni_agd_arm + ni_agd_contrast))] : Xdummy;

// Split offsets into IPD and AgD rows
vector[0] odummy;
vector[has_offset && ni_ipd ? ni_ipd : 0] offset_ipd = has_offset && ni_ipd ? offsets[1:ni_ipd] : odummy;
vector[has_offset && ni_agd_arm ? nint_max * ni_agd_arm : 0] offset_agd_arm = has_offset && ni_agd_arm ? offsets[(ni_ipd + 1):(ni_ipd + nint_max * ni_agd_arm)] : odummy;
vector[has_offset && ni_agd_contrast ? nint_max * ni_agd_contrast : 0] offset_agd_contrast = has_offset && ni_agd_contrast ? offsets[(ni_ipd + nint_max * ni_agd_arm + 1):(ni_ipd + nint_max * (ni_agd_arm + ni_agd_contrast))] : odummy;

// nint/int_thin for numerical integration checks
int n_int_thin = (nint_max > 1 && int_thin > 0) ? nint %/% int_thin : 0;

// Inverse covariance matrix for contrasts
matrix[ni_agd_contrast ? ni_agd_contrast : 1, ni_agd_contrast ? ni_agd_contrast : 1] inv_Sigma = inverse_spd(agd_contrast_Sigma);

// Construct number of contrasts in each study for contrast-based AgD by looking at Sigma covariance matrix
// NOTE: Sigma must be block diagonal (i.e. all contrasts for a single study together)
array[ns_agd_contrast] int nc_agd_contrast;
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

// Number of class effects and which_class_sd, mapping class sd to classes
int<lower=0> n_class = class_effects ? max(which_CE) : 0;
array[class_effects ? nt-1 : 0] int<lower=0> which_class; // Vector to mapping treatments to class effects
array[class_effects ? n_class : 0] int<lower=0> which_class_sd; // Vector to map class SDs to classes

int<lower=0> n_class_trts = num_elements(which_gt0a(which_CE));  // Number of treatments with class effects
array[n_class_trts] int which_class_trt = which_gt0a(which_CE); // Vector mapping classes to treatments
 if (class_effects) {
    for (c in 1:n_class) {
      for (t in 1:nt - 1) {
        if (which_CE[t] == c) {
          which_class_sd[c] = which_CE_sd[t];
          break;
        }
      }
    }
  }

if (class_effects) {
  int count = 1;  // Counter for nonzero elements in f_class

  for (t in 1:(nt - 1)) {
    if (which_CE_sd[t] > 0) {
      which_class[t] = count;
      count += 1;
    } else {
      which_class[t] = 0;  // No class effect for this treatment
    }
  }
}
