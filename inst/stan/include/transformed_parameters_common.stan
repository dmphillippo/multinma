// Common definitions for the transformed parameters block

// -- Likelihood parameters needed later for log lik calculation --
vector[ni_ipd] eta_ipd; // IPD linear predictor
vector[ni_ipd] theta_ipd; // IPD transformed predictor

// -- RE deltas --
vector[n_delta] f_delta = sigma * RE_L * u_delta;

// -- Back-transformed parameters --
vector[nX] allbeta = QR ? R_inv * beta_tilde : beta_tilde;
// Study baselines
vector[ns_ipd + ns_agd_arm + ns_agd_contrast] beta0 = allbeta[1:totns];
// Treatment effects
vector[nt - 1] gamma = allbeta[(totns +1):(totns + nt - 1)];
// Regression predictors
vector[nPV] beta = allbeta[(totns + nt):];

// -- AgD integration --
vector[nint * (ni_agd_arm + ni_agd_contrast)] theta_ii;
vector[ni_agd_arm + ni_agd_contrast] theta_bar;

// -- IPD model --
// We define the IPD and AgD models here in the transformed parameters block,
// as the linear predictors are required to calculate the log likelihood
// later on. This is slightly more inefficient than defining the models
// locally in the model block.
{
  vector[ni_ipd] eta_ipd_noRE = Q_ipd * beta_tilde;
  for (i in 1:ni_ipd) {
    if (delta_design[ipd_arm[i]])
      eta_ipd[i] = eta_ipd_noRE[i] + f_delta[delta_design[ipd_arm[i]]];
    else
      eta_ipd[i] = eta_ipd_noRE[i];
  }
}
