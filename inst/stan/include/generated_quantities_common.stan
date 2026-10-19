// Common definitions for the generated quantities block

// -- Log likelihood, residual deviance, fitted values --
vector[ni_ipd + ni_agd_arm + ns_agd_contrast + ns_agd_regression] log_lik;
vector[ni_ipd + ni_agd_arm + ns_agd_contrast + ns_agd_regression] resdev;
vector[ni_agd_contrast] fitted_agd_contrast;
vector[nc_agd_regression] fitted_agd_regression;

// -- Estimate integration error --
  vector[ni_agd_contrast * n_int_thin] theta_bar_cum_agd_contrast;

// -- RE shrunken estimate delta --
// Note: These are the individual-level trial-specific treatment effects
vector[n_delta] delta;

if (RE) {
  int s = 1;
  for (i in 1:(narm_ipd + narm_agd_arm)) {
    if (which_RE[i]) {
      delta[s] = (trt[i] > 1 ? d[trt[i] - 1] : 0) + (which_RE[i] ? f_delta[which_RE[i]] : 0);
      s += 1;
    }
  }
  for (i in 1:ni_agd_contrast) {
    if (which_RE[narm_ipd + narm_agd_arm + i]) {
      delta[s] = (trt[narm_ipd + narm_agd_arm + i] > 1 ? d[trt[narm_ipd + narm_agd_arm + i] - 1] : 0) +
        (which_RE[narm_ipd + narm_agd_arm + i] ? f_delta[which_RE[narm_ipd + narm_agd_arm + i]] : 0);

      if (agd_contrast_trt_b[i] > 1)
        delta[s] -= d[agd_contrast_trt_b[i] - 1];

      s += 1;
    }
  }
}

// Integration error for AgD
for (i in 1:ni_agd_contrast) {
  for (j in 1:n_int_thin) {
    theta_bar_cum_agd_contrast[(i-1)*n_int_thin + j] = mean(eta_agd_contrast_ii[(1 + (i-1)*nint):((i-1)*nint + j*int_thin)]);
  }
}
// Fitted values
fitted_agd_contrast = eta_agd_contrast_bar;
fitted_agd_regression = eta_agd_regression;

// log likelihood, residual deviance for contrast-based AgD are *per study*
{
  int a = 0;
  int nc;
  for (s in 1:ns_agd_contrast) {
    nc = nc_agd_contrast[s];
    log_lik[ni_ipd + ni_agd_arm + s] = multi_normal_lpdf(agd_contrast_y[(a + 1):(a + nc)] |
      eta_agd_contrast_bar[(a + 1):(a + nc)], agd_contrast_Sigma[(a + 1):(a + nc), (a + 1):(a + nc)]);
    resdev[ni_ipd + ni_agd_arm + s] = quad_form(inv_Sigma[(a + 1):(a + nc), (a + 1):(a + nc)],
                                                agd_contrast_y[(a + 1):(a + nc)] - eta_agd_contrast_bar[(a + 1):(a + nc)]);
    a += nc;
  }
}

// log likelihood, residual deviance for regression-based AgD are *per study*
{
  int c=0;
  for (i in 1:ns_agd_regression) {
    vector[agd_regression_ncoef[i]] est = segment(agd_regression_est, c + 1, agd_regression_ncoef[i]);
    vector[agd_regression_ncoef[i]] eta = segment(eta_agd_regression, c + 1, agd_regression_ncoef[i]);
    matrix[agd_regression_ncoef[i], agd_regression_ncoef[i]] sigma_chol = block(agd_regression_cov[i], 1, 1, agd_regression_ncoef[i], agd_regression_ncoef[i]);
    matrix[agd_regression_ncoef[i], agd_regression_ncoef[i]] sigma_inv = chol2inv(sigma_chol);

    log_lik[ni_ipd + ni_agd_arm + ns_agd_contrast + i] = multi_normal_cholesky_lpdf(est | eta, sigma_chol);
    resdev[ni_ipd + ni_agd_arm + ns_agd_contrast + i] = quad_form(sigma_inv, est - eta);
    c += agd_regression_ncoef[i];
  }
}

