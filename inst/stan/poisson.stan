functions {
#include /include/prior_select.stan
#include /include/count_nonzero.stan
}
data {
#include /include/data_common.stan

  // Outcomes
  int<lower=0> ipd_r[ni_ipd];
  vector<lower=0>[ni_ipd] ipd_E;
  int<lower=0> agd_arm_r[ni_agd_arm];
  vector<lower=0>[ni_agd_arm] agd_arm_E;
}
transformed data {
#include /include/transformed_data_common.stan
}
parameters {
#include /include/parameters_common.stan
}
transformed parameters {
  vector<lower=0>[ni_ipd] E_eta_ipd;
  vector<lower=0>[ni_agd_arm] E_theta_agd_arm;

#include /include/transformed_parameters_theta.stan
#include /include/transformed_parameters_common.stan

  // -- IPD model --
  if (link == 1) // log link
    theta_ipd = exp(eta_ipd);

  // -- AgD model (arm-based) --
  if (ni_agd_arm) {
    if (nint > 1) { // -- If integration points are used --
      if (RE) {
        vector[nint * ni_agd_arm] eta_agd_arm_noRE = has_offset ?
          X_agd_arm * beta_tilde + offset_agd_arm :
          X_agd_arm * beta_tilde;

        if (link == 1) { // log link
          for (i in 1:ni_agd_arm) {
            if (which_RE[narm_ipd + i])
              theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint)] = exp(eta_agd_arm_noRE[(1 + (i-1)*nint):(i*nint)] + f_delta[which_RE[narm_ipd + i]]);
            else
              theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint)] = exp(eta_agd_arm_noRE[(1 + (i-1)*nint):(i*nint)]);
          }
        }

        for (i in 1:ni_agd_arm) {
          theta_agd_arm_bar[i] = mean(theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint)]);
        }

      } else {
        if (link == 1) { // log link
          theta_agd_arm_ii = has_offset ?
            exp(X_agd_arm * beta_tilde + offset_agd_arm) :
            exp(X_agd_arm * beta_tilde);
        }

        for (i in 1:ni_agd_arm) {
          theta_agd_arm_bar[i] = mean(theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint)]);
        }
      }
    } else { // -- If no integration --
      if (RE) {
        vector[nint * ni_agd_arm] eta_agd_arm_noRE = has_offset ?
          X_agd_arm * beta_tilde + offset_agd_arm :
          X_agd_arm * beta_tilde;

        if (link == 1) { // log link
          for (i in 1:ni_agd_arm) {
            if (which_RE[narm_ipd + i])
              theta_agd_arm_bar[i] = exp(eta_agd_arm_noRE[i] + f_delta[which_RE[narm_ipd + i]]);
            else
              theta_agd_arm_bar[i] = exp(eta_agd_arm_noRE[i]);
          }
        }
      } else {
        if (link == 1) // log link
          theta_agd_arm_bar = has_offset ?
            exp(X_agd_arm * beta_tilde + offset_agd_arm) :
            exp(X_agd_arm * beta_tilde);
      }
    }
  }

  // Predictors with time at risk offset
  E_eta_ipd = eta_ipd + ipd_E;
  E_theta_agd_arm = theta_agd_arm_bar .* agd_arm_E;
}
model {
#include /include/model_common.stan

  // -- IPD likelihood --
  // Could replace log link sampling statement with poisson_log_glm in Stan > 2.20
  if (link == 1) { // log link
    ipd_r ~ poisson_log(E_eta_ipd);
  }

  // -- AgD likelihood (arm-based) --
  agd_arm_r ~ poisson(E_theta_agd_arm);
}
generated quantities {
#include /include/generated_quantities_theta_fitted.stan
#include /include/generated_quantities_common.stan
#include /include/generated_quantities_theta.stan

  // IPD log likelihood and residual deviance
  if (link == 1) { // log link
    vector[ni_ipd] E_theta_ipd = exp(E_eta_ipd);
    for (i in 1:ni_ipd) {
      log_lik[i] = poisson_log_lpmf(ipd_r[i] | E_eta_ipd[i]);
      resdev[i] = 2 * ((E_theta_ipd[i] - ipd_r[i]) + lmultiply(ipd_r[i], ipd_r[i] / E_theta_ipd[i]));
      fitted_ipd[i] = E_theta_ipd[i];
    }
  }

  // AgD (arm-based) log likelihood and residual deviance
  for (i in 1:ni_agd_arm) {
    log_lik[ni_ipd + i] = poisson_lpmf(agd_arm_r[i] | E_theta_agd_arm[i]);
    resdev[ni_ipd + i] = 2 * ((E_theta_agd_arm[i] - agd_arm_r[i]) +
                              lmultiply(agd_arm_r[i], agd_arm_r[i] / E_theta_agd_arm[i]));
    fitted_agd_arm[i] = E_theta_agd_arm[i];
  }

}
