functions {
#include /include/prior_select.stan
}
data {
#include /include/data_common.stan

  // Outcomes
  int<lower=0, upper=1> ipd_r[ni_ipd];
  int<lower=0> agd_arm_n[ni_agd_arm];
  int<lower=0> agd_arm_r[ni_agd_arm];
}
transformed data {
#include /include/transformed_data_common.stan
}
parameters {
#include /include/parameters_common.stan
}
transformed parameters {
#include /include/transformed_parameters_common.stan

  // -- IPD model --
  if (link == 1) // logit link
    theta_ipd = inv_logit(eta_ipd);
  else if (link == 2) // probit link
    theta_ipd = Phi(eta_ipd);

  // -- AgD model (arm-based) --
  if (ni_agd_arm) {
    if (nint > 1) { // -- If integration points are used --
      if (RE) {
        vector[nint * ni_agd_arm] eta_agd_arm_noRE = X_agd_arm * beta_tilde;

        if (link == 1) { // logit link
          for (i in 1:ni_agd_arm) {
            if (which_RE[narm_ipd + i])
              theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint)] = inv_logit(eta_agd_arm_noRE[(1 + (i-1)*nint):(i*nint)] + f_delta[which_RE[narm_ipd + i]]);
            else
              theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint)] = inv_logit(eta_agd_arm_noRE[(1 + (i-1)*nint):(i*nint)]);
          }
        } else if (link == 2) { // probit link
          for (i in 1:ni_agd_arm) {
            if (which_RE[narm_ipd + i])
              theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint)] = Phi(eta_agd_arm_noRE[(1 + (i-1)*nint):(i*nint)] + f_delta[which_RE[narm_ipd + i]]);
            else
              theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint)] = Phi(eta_agd_arm_noRE[(1 + (i-1)*nint):(i*nint)]);
          }
        }

        for (i in 1:ni_agd_arm) {
          theta_agd_arm_bar[i] = mean(theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint)]);
        }

      } else {
        if (link == 1) { // logit link
          theta_agd_arm_ii = inv_logit(X_agd_arm * beta_tilde);
        } else if (link == 2) { // probit link
          theta_agd_arm_ii = Phi(X_agd_arm * beta_tilde);
        }

        for (i in 1:ni_agd_arm) {
          theta_agd_arm_bar[i] = mean(theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint)]);
        }
      }
    } else { // -- If no integration --
      if (RE) {
        vector[nint * ni_agd_arm] eta_agd_arm_noRE = X_agd_arm * beta_tilde;

        if (link == 1) { // logit link
          for (i in 1:ni_agd_arm) {
            if (which_RE[narm_ipd + i])
              theta_agd_arm_bar[i] = inv_logit(eta_agd_arm_noRE[i] + f_delta[which_RE[narm_ipd + i]]);
            else
              theta_agd_arm_bar[i] = inv_logit(eta_agd_arm_noRE[i]);
          }
        } else if (link == 2) { // probit link
          for (i in 1:ni_agd_arm) {
            if (which_RE[narm_ipd + i])
              theta_agd_arm_bar[i] = Phi(eta_agd_arm_noRE[i] + f_delta[which_RE[narm_ipd + i]]);
            else
              theta_agd_arm_bar[i] = Phi(eta_agd_arm_noRE[i]);
          }
        }
      } else {
        if (link == 1) // logit link
          theta_agd_arm_bar = inv_logit(X_agd_arm * beta_tilde);
        else if (link == 2) // probit link
          theta_agd_arm_bar = Phi(X_agd_arm * beta_tilde);
      }
    }
  }
}
model {
#include /include/model_common.stan

  // -- IPD likelihood --
  if (link == 1) { // logit link
    // Could replace with bernoulli_logit_glm in Stan > 2.20
    ipd_r ~ bernoulli_logit(eta_ipd);
  } else {
    ipd_r ~ bernoulli(theta_ipd);
  }

  // -- AgD likelihood (arm-based) --
  agd_arm_r ~ binomial(agd_arm_n, theta_agd_arm_bar);
}
generated quantities {
#include /include/generated_quantities_common.stan

  // IPD log likelihood and residual deviance
  for (i in 1:ni_ipd) {
    log_lik[i] = bernoulli_lpmf(ipd_r[i] | theta_ipd[i]);
    resdev[i] = -2 * log_lik[i];
  }

  // AgD (arm-based) log likelihood and residual deviance
  for (i in 1:ni_agd_arm) {
    log_lik[ni_ipd + i] = binomial_lpmf(agd_arm_r[i] | agd_arm_n[i], theta_agd_arm_bar[i]);
    resdev[ni_ipd + i] = 2 *
      (lmultiply(agd_arm_r[i],
                 agd_arm_r[i] / (agd_arm_n[i] * theta_agd_arm_bar[i])) +
       lmultiply(agd_arm_n[i] - agd_arm_r[i],
                 (agd_arm_n[i] - agd_arm_r[i]) / (agd_arm_n[i] - agd_arm_n[i] * theta_agd_arm_bar[i])));

	  for (j in 1:n_int_thin) {
      theta_bar_cum[(i-1)*n_int_thin + j] = mean(theta_agd_arm_ii[(1 + (i-1)*nint):((i-1)*nint + j*int_thin)]);
    }
  }

}
