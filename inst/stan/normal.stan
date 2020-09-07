functions {
#include /include/prior_select.stan
#include /include/count_nonzero.stan
}
data {
#include /include/data_common.stan

  // Prior on IPD variance
  int<lower=0,upper=5> prior_aux_dist;
  real prior_aux_location;
  real<lower=0> prior_aux_scale;
  real<lower=0> prior_aux_df;

  // Outcomes
  real ipd_y[ni_ipd];
  real agd_arm_y[ni_agd_arm];
  real<lower=0> agd_arm_se[ni_agd_arm];
}
transformed data {
#include /include/transformed_data_common.stan
}
parameters {
#include /include/parameters_common.stan
  vector<lower=0>[narm_ipd] sigma;
}
transformed parameters {
#include /include/transformed_parameters_theta.stan
#include /include/transformed_parameters_common.stan

  // -- IPD model --
  if (link == 1) // identity link
    theta_ipd = eta_ipd;
  else if (link == 2) // log link
    theta_ipd = exp(eta_ipd);

  // -- AgD model (arm-based) --
  if (ni_agd_arm) {
    if (nint > 1) { // -- If integration points are used --
      if (RE) {
        vector[nint * ni_agd_arm] eta_agd_arm_noRE = has_offset ?
          X_agd_arm * beta_tilde + offset_agd_arm :
          X_agd_arm * beta_tilde;

        if (link == 1) { // identity link
          for (i in 1:ni_agd_arm) {
            if (which_RE[narm_ipd + i])
              theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint)] = eta_agd_arm_noRE[(1 + (i-1)*nint):(i*nint)] + f_delta[which_RE[narm_ipd + i]];
            else
              theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint)] = eta_agd_arm_noRE[(1 + (i-1)*nint):(i*nint)];
          }
        } else if (link == 2) { // log link
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
        if (link == 1) { // identity link
          theta_agd_arm_ii = has_offset ?
            X_agd_arm * beta_tilde + offset_agd_arm :
            X_agd_arm * beta_tilde;
        } else if (link == 2) { // log link
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

        if (link == 1) { // identity link
          for (i in 1:ni_agd_arm) {
            if (which_RE[narm_ipd + i])
              theta_agd_arm_bar[i] = eta_agd_arm_noRE[i] + f_delta[which_RE[narm_ipd + i]];
            else
              theta_agd_arm_bar[i] = eta_agd_arm_noRE[i];
          }
        } else if (link == 2) { // log link
          for (i in 1:ni_agd_arm) {
            if (which_RE[narm_ipd + i])
              theta_agd_arm_bar[i] = exp(eta_agd_arm_noRE[i] + f_delta[which_RE[narm_ipd + i]]);
            else
              theta_agd_arm_bar[i] = exp(eta_agd_arm_noRE[i]);
          }
        }
      } else {
        if (link == 1) // identity link
          theta_agd_arm_bar = has_offset ?
            X_agd_arm * beta_tilde + offset_agd_arm :
            X_agd_arm * beta_tilde;
        else if (link == 2) // log link
          theta_agd_arm_bar = has_offset ?
            exp(X_agd_arm * beta_tilde + offset_agd_arm) :
            exp(X_agd_arm * beta_tilde);
      }
    }
  }
}
model {
#include /include/model_common.stan

  // -- Prior on arm-level variance --
  prior_select_lp(sigma, prior_aux_dist, prior_aux_location, prior_aux_scale, prior_aux_df);

  // -- IPD likelihood --
  // Could replace identity link sampling statement with normal_id_glm in Stan > 2.20
  ipd_y ~ normal(theta_ipd, sigma[ipd_arm]);

  // -- AgD likelihood (arm-based) --
  agd_arm_y ~ normal(theta_agd_arm_bar, agd_arm_se);
}
generated quantities {
#include /include/generated_quantities_theta_fitted.stan
#include /include/generated_quantities_common.stan
#include /include/generated_quantities_theta.stan

  // IPD log likelihood and residual deviance
  for (i in 1:ni_ipd) {
    log_lik[i] = normal_lpdf(ipd_y[i] | theta_ipd[i], sigma[ipd_arm[i]]);
    fitted_ipd[i] = theta_ipd[i];
    resdev[i] = (ipd_y[i] - fitted_ipd[i])^2 / sigma[ipd_arm[i]]^2;
  }

  // AgD (arm-based) log likelihood and residual deviance
  for (i in 1:ni_agd_arm) {
    log_lik[ni_ipd + i] = normal_lpdf(agd_arm_y[i] | theta_agd_arm_bar[i], agd_arm_se[i]);
    fitted_agd_arm[i] = theta_agd_arm_bar[i];
    resdev[ni_ipd + i] = (agd_arm_y[i] - fitted_agd_arm[i])^2 / agd_arm_se[i]^2;
  }

}
