functions {
  #include prior_select.stan
}
data {
  #include data_common.stan

  // Outcomes
  int<lower=0, upper=1> ipd_r[ni_ipd];
  int<lower=0> agd_arm_n[ni_agd];
  int<lower=0> agd_arm_r[ni_agd];
}
transformed data {
  #include transformed_data_common.stan
}
parameters {
  #include parameters_common.stan
  real<lower=0> nprime[ni_agd_arm];
  real<lower=0, upper=1> pprime[ni_agd_arm];
}
transformed parameters {
  #include transformed_parameters_common.stan

  // -- IPD model --
  if (link == 1) // logit link
    theta_ipd = inv_logit(eta_ipd);
  else if (link == 2) // probit link
    theta_ipd = Phi(eta_ipd);

  // -- AgD model (arm-based) --
  {
    vector[nint * ni_agd_arm] eta_agd_arm_noRE = Q_agd_arm * beta_tilde;

    if (link == 1) { // logit link
      for (i in 1:ni_agd_arm) {
        if (delta_design[narm_ipd + i])
          theta_ii[(1 + (i-1)*nint):(i*nint)] = inv_logit(eta_agd_noRE[(1 + (i-1)*nint):(i*nint)] + f_delta[delta_design[narm_ipd + i]]);
        else
          theta_ii[(1 + (i-1)*nint):(i*nint)] = inv_logit(eta_agd_noRE[(1 + (i-1)*nint):(i*nint)]);
      }
    } else if (link == 2) { // probit link
      for (i in 1:ni_agd_arm) {
        if (delta_design[narm_ipd + i])
          theta_ii[(1 + (i-1)*nint):(i*nint)] = Phi(eta_agd_noRE[(1 + (i-1)*nint):(i*nint)] + f_delta[delta_design[narm_ipd + i]]);
        else
          theta_ii[(1 + (i-1)*nint):(i*nint)] = Phi(eta_agd_noRE[(1 + (i-1)*nint):(i*nint)]);
      }
    }

    for (i in 1:ni_agd_arm) {
      theta_bar[i] = mean(theta_ii[(1 + (i-1)*nint):(i*nint)]);
      theta2_bar[i] = dot_self(theta_ii[(1 + (i-1)*nint):(i*nint)]) / nint;

      // Calculate adjusted n and p
      nprime[i] = ag_n[i] * theta_bar[i]^2 / theta2_bar[i];
      pprime[i] = theta2_bar[i] / theta_bar[i];

      // Reject sample if nprime less than number of observed events (shouldn't be necessary...)
      if (nprime[i] < agd_arm_r[i]) reject("nprime = ", nprime[i], " less than r = ", agd_arm_r[i]);
    }
  }
}
model {
  #include priors.stan

  // -- Random effects --
  u_delta ~ normal(0, 1);

  // -- IPD likelihood --
  if (link == 1) { // logit link
    // Could replace with bernoulli_logit_glm in Stan > 2.20
    y ~ bernoulli_logit(eta_ipd);
  } else {
    y ~ bernoulli(theta_ipd);
  }

  // -- AgD likelihood (arm-based) --
  // We have to hand code the log likelihood contribution for the adjusted
  // binomial here, as N is not necessarily an integer (which Stan doesn't
  // like). The following is exactly equivalent to:
  // ag_r ~ binomial(nprime, pprime);
  for (i in 1:ni_agd)
    target += lchoose(nprime[i], agd_arm_r[i]) +
              lmultiply(agd_arm_r[i], pprime[i]) +
              (nprime[i] - agd_arm_r[i]) * log1m(pprime[i]);
}
generated quantities {
  // -- Log likelihood and residual deviance calculation --
  vector[ni_ipd + ni_agd_arm + ni_agd_contrast] log_lik;
  vector[ni_ipd + ni_agd_arm + ni_agd_contrast] resdev;

  // -- Estimate integration error --
  vector[ni_agd * n_int_thin] theta_bar_cum;
  vector[ni_agd * n_int_thin] theta2_bar_cum;

  // -- RE shrunken estimate delta --
  // Note: These are the individual-level trial-specific treatment effects
  vector[narm_ipd + ni_agd_arm + ni_agd_contrast] delta;

  // For the shrunken estimates, since trt 1 is the reference and REs are treatment based
  // rather than arm based, any trt 1 arm has delta = 0
  for (i in 1:(narm_ipd + ni_agd_arm + ni_agd_contrast)) {
    delta[i] = trt[i] == 1 ? 0 : gamma[trt[i] - 1] + f_delta[delta_design[i]];
  }

  for (i in 1:ni_ipd) {
	  p_hat[i] = theta[i];
	  r_hat[i] = theta[i];
    log_lik[i] = bernoulli_lpmf(y[i] | theta[i]);
    resdev[i] = -2 * log_lik[i];
  }

  for (i in 1:ni_agd) {
    log_lik[ni_ipd + i] = lchoose(nprime[i], ag_r[i]) + lmultiply(ag_r[i], pprime[i]) + (nprime[i] - ag_r[i]) * log1m(pprime[i]);
    // p_bar_diff[i] = p_bar[i] - mean(p_ii[(1 + (i-1)*nint):((i-1)*nint + half_nint)]);
    // p2_bar_diff[i] = p2_bar[i] - (dot_self(p_ii[(1 + (i-1)*nint):((i-1)*nint + half_nint)]) / (half_nint));

    r_hat[ni_ipd + i] = nprime[i] * pprime[i];
    p_hat[ni_ipd + i] = r_hat[i] / ag_n[i];

    // Approximate residual deviance for AgD, letting nprime be fixed
    resdev[ni_ipd + i] = 2 * (lmultiply(ag_r[i], ag_r[i] / (nprime[i] * pprime[i])) + lmultiply(ag_n[i] - ag_r[i], (ag_n[i] - ag_r[i]) / (ag_n[i] - nprime[i] * pprime[i])));

	for (j in 1:n_int_thin) {
      p_bar_cum[(i-1)*n_int_thin + j] = mean(p_ii[(1 + (i-1)*nint):((i-1)*nint + j*int_thin)]);
      p2_bar_cum[(i-1)*n_int_thin + j] = (dot_self(p_ii[(1 + (i-1)*nint):((i-1)*nint + j*int_thin)]) / (j*int_thin));
    }
  }

}
