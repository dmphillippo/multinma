// Common definitions for the transformed parameters block

// -- Linear and transformed predictors --
vector[ni_ipd] eta_ipd; // IPD linear predictor
vector[ni_ipd] theta_ipd; // IPD transformed predictor

vector[ni_agd_arm] theta_agd_arm; // AgD (arm-based) transformed predictor

vector[ni_agd_contrast] eta_agd_contrast; // AgD (contrast-based) linear predictor

// -- RE deltas --
vector[n_delta] f_delta = tau[1] * RE_L * u_delta;

// -- Back-transformed parameters --
vector[nX] allbeta = QR ? R_inv * beta_tilde : beta_tilde;
// Study baselines
vector[totns] mu = allbeta[1:totns];
// Treatment effects
vector[nt - 1] gamma = allbeta[(totns +1):(totns + nt - 1)];
// Regression predictors
vector[nX] beta = allbeta[(totns + nt):];

// -- AgD integration --
vector[nint * ni_agd_arm] theta_agd_arm_ii;
vector[ni_agd_arm] theta_agd_arm_bar;
vector[nint * ni_agd_contrast] eta_agd_contrast_ii;
vector[ni_agd_contrast] eta_agd_contrast_bar;

// -- IPD model --
// We define the IPD and AgD models here in the transformed parameters block,
// as the linear predictors are required to calculate the log likelihood
// later on. This is slightly more inefficient than defining the models
// locally in the model block.
{
  vector[ni_ipd] eta_ipd_noRE = X_ipd * beta_tilde;
  for (i in 1:ni_ipd) {
    if (which_RE[ipd_arm[i]])
      eta_ipd[i] = eta_ipd_noRE[i] + f_delta[which_RE[ipd_arm[i]]];
    else
      eta_ipd[i] = eta_ipd_noRE[i];
  }
}

// -- AgD model (contrast-based) --
{
  vector[nint * ni_agd_contrast] eta_agd_contrast_noRE = X_agd_contrast * beta_tilde;
  for (i in 1:ni_agd_contrast) {
    if (which_RE[narm_ipd + ni_agd_arm + i])
      eta_agd_contrast_ii[(1 + (i-1)*nint):(i*nint)] =
        eta_agd_contrast_noRE[(1 + (i-1)*nint):(i*nint)] + f_delta[which_RE[narm_ipd + ni_agd_arm + i]];
    else
      eta_agd_contrast_ii[(1 + (i-1)*nint):(i*nint)] =
        eta_agd_contrast_noRE[(1 + (i-1)*nint):(i*nint)];

    eta_agd_contrast_bar[i] = mean(eta_agd_contrast_ii[(1 + (i-1)*nint):(i*nint)]);
  }
}
