functions {
#include /include/prior_select.stan
#include /include/count_nonzero.stan
}
data {
#include /include/data_common.stan

  // Outcomes
  array[ni_ipd] int<lower=0, upper=1> ipd_r;
  array[ni_agd_arm] int<lower=0> agd_arm_n;
  array[ni_agd_arm] int<lower=0> agd_arm_r;
}
transformed data {
#include /include/transformed_data_common.stan
}
parameters {
#include /include/parameters_common.stan
}
transformed parameters {
  vector[ni_agd_arm] theta2_agd_arm_bar;
  array[ni_agd_arm] real<lower=0> nprime;
  array[ni_agd_arm] real<lower=0, upper=1> pprime;
#include /include/transformed_parameters_theta.stan
#include /include/transformed_parameters_common.stan

  // -- IPD model --
  if (link == 1) // logit link
    theta_ipd = inv_logit(eta_ipd);
  else if (link == 2) // probit link
    theta_ipd = Phi(eta_ipd);
  else if (link == 3) // cloglog link
    theta_ipd = inv_cloglog(eta_ipd);

  // -- AgD model (arm-based) --
  if (ni_agd_arm) {
    if (nint_max > 1) { // -- If integration points are used --
      vector[nint_max * ni_agd_arm] eta_agd_arm_noRE = has_offset ?
          X_agd_arm * beta_tilde + offset_agd_arm :
          X_agd_arm * beta_tilde;

    if (class_effects) {
      for (i in 1:ni_agd_arm) {
        if (agd_arm_trt[i] > 1 && which_CE[agd_arm_trt[i] - 1]) {
          eta_agd_arm_noRE[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)] += f_class[which_class[agd_arm_trt[i] - 1]];
        }
      }
    }

      if (RE) {
        if (link == 1) { // logit link
          for (i in 1:ni_agd_arm) {
            if (which_RE[narm_ipd + i])
              theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint)] = inv_logit(eta_agd_arm_noRE[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)] + f_delta[which_RE[narm_ipd + i]]);
            else
              theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint)] = inv_logit(eta_agd_arm_noRE[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)]);
          }
        } else if (link == 2) { // probit link
          for (i in 1:ni_agd_arm) {
            if (which_RE[narm_ipd + i])
              theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint)] = Phi(eta_agd_arm_noRE[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)] + f_delta[which_RE[narm_ipd + i]]);
            else
              theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint)] = Phi(eta_agd_arm_noRE[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)]);
          }
        } else if (link == 3) { // cloglog link
          for (i in 1:ni_agd_arm) {
            if (which_RE[narm_ipd + i])
              theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint)] = inv_cloglog(eta_agd_arm_noRE[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)] + f_delta[which_RE[narm_ipd + i]]);
            else
              theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint)] = inv_cloglog(eta_agd_arm_noRE[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)]);
          }
        }
      } else {
        if (link == 1) { // logit link
          if (nint == nint_max) theta_agd_arm_ii = inv_logit(eta_agd_arm_noRE);
          else for (i in 1:ni_agd_arm) theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint)] = inv_logit(eta_agd_arm_noRE[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)]);
        } else if (link == 2) { // probit link
          if (nint == nint_max) theta_agd_arm_ii = Phi(eta_agd_arm_noRE);
          else for (i in 1:ni_agd_arm) theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint)] = Phi(eta_agd_arm_noRE[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)]);
        } else if (link == 3) { // cloglog link
          if (nint == nint_max) theta_agd_arm_ii = inv_cloglog(eta_agd_arm_noRE);
          else for (i in 1:ni_agd_arm) theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint)] = inv_cloglog(eta_agd_arm_noRE[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)]);
        }
      }

      for (i in 1:ni_agd_arm) {
        theta_agd_arm_bar[i] = mean(theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint)]);
        theta2_agd_arm_bar[i] = dot_self(theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint)]) / nint;

        // Calculate adjusted n and p
        nprime[i] = agd_arm_n[i] * theta_agd_arm_bar[i]^2 / theta2_agd_arm_bar[i];
        pprime[i] = theta2_agd_arm_bar[i] / theta_agd_arm_bar[i];

        // Reject sample if nprime less than number of observed events (shouldn't be necessary...)
        if (nprime[i] < agd_arm_r[i]) reject("nprime = ", nprime[i], " less than r = ", agd_arm_r[i]);
      }
    } else { // -- If no integration --
      if (RE) {
        vector[nint * ni_agd_arm] eta_agd_arm_noRE = has_offset ?
          X_agd_arm * beta_tilde + offset_agd_arm :
          X_agd_arm * beta_tilde;

      // Add class effects contribution to the linear predictor
        if (class_effects) {
          for (i in 1:ni_agd_arm) {
            if (agd_arm_trt[i] > 1 && which_CE[agd_arm_trt[i] - 1]) {
              eta_agd_arm_noRE[i] += f_class[which_class[agd_arm_trt[i] - 1]];
            }
          }
        }

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
        } else if (link == 3) { // cloglog link
          for (i in 1:ni_agd_arm) {
            if (which_RE[narm_ipd + i])
              theta_agd_arm_bar[i] = inv_cloglog(eta_agd_arm_noRE[i] + f_delta[which_RE[narm_ipd + i]]);
            else
              theta_agd_arm_bar[i] = inv_cloglog(eta_agd_arm_noRE[i]);
          }
        }
      } else {

        vector[ni_agd_arm] eta_agd_arm_noRE = has_offset ?
          X_agd_arm * beta_tilde + offset_agd_arm :
          X_agd_arm * beta_tilde;

        // Add class effects contribution to the linear predictor
        if (class_effects) {
          for (i in 1:ni_agd_arm) {
            if (agd_arm_trt[i] > 1 && which_CE[agd_arm_trt[i] - 1]) {
              eta_agd_arm_noRE[i] += f_class[which_class[agd_arm_trt[i] - 1]];
            }
          }
        }

        if (link == 1) // logit link
          theta_agd_arm_bar = inv_logit(eta_agd_arm_noRE) ;
        else if (link == 2) // probit link
          theta_agd_arm_bar = Phi(eta_agd_arm_noRE) ;
        else if (link == 3) // cloglog link
          theta_agd_arm_bar = inv_cloglog(eta_agd_arm_noRE) ;
      }

      theta2_agd_arm_bar = theta_agd_arm_bar .* theta_agd_arm_bar;

      for (i in 1:ni_agd_arm) {
        // Calculate adjusted n and p
        nprime[i] = agd_arm_n[i] * theta_agd_arm_bar[i]^2 / theta2_agd_arm_bar[i];
        pprime[i] = theta2_agd_arm_bar[i] / theta_agd_arm_bar[i];

        // Reject sample if nprime less than number of observed events (shouldn't be necessary...)
        if (nprime[i] < agd_arm_r[i]) reject("nprime = ", nprime[i], " less than r = ", agd_arm_r[i]);
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
  // We have to hand code the log likelihood contribution for the adjusted
  // binomial here, as N is not necessarily an integer (which Stan doesn't
  // like). The following is exactly equivalent to:
  // ag_r ~ binomial(nprime, pprime);
  for (i in 1:ni_agd_arm)
    target += lchoose(nprime[i], agd_arm_r[i]) +
              lmultiply(agd_arm_r[i], pprime[i]) +
              (nprime[i] - agd_arm_r[i]) * log1m(pprime[i]);
}
generated quantities {
  vector[ni_agd_arm * n_int_thin] theta2_bar_cum;
#include /include/generated_quantities_theta_fitted.stan
#include /include/generated_quantities_common.stan
#include /include/generated_quantities_theta.stan

  // IPD log likelihood and residual deviance
  for (i in 1:ni_ipd) {
    log_lik[i] = bernoulli_lpmf(ipd_r[i] | theta_ipd[i]);
    resdev[i] = -2 * log_lik[i];
    fitted_ipd[i] = theta_ipd[i];
  }

  // AgD (arm-based) log likelihood and residual deviance
  for (i in 1:ni_agd_arm) {
    log_lik[ni_ipd + i] = lchoose(nprime[i], agd_arm_r[i]) +
                          lmultiply(agd_arm_r[i], pprime[i]) +
                          (nprime[i] - agd_arm_r[i]) * log1m(pprime[i]);
    // Approximate residual deviance for AgD, letting nprime be fixed
    resdev[ni_ipd + i] = 2 * ((agd_arm_r[i] > 0 ?
                                lmultiply(agd_arm_r[i],
                                          agd_arm_r[i] / (nprime[i] * pprime[i])) : 0) +
                              (agd_arm_r[i] < agd_arm_n[i] ?
                                lmultiply(agd_arm_n[i] - agd_arm_r[i],
                                          (agd_arm_n[i] - agd_arm_r[i]) / (agd_arm_n[i] - nprime[i] * pprime[i])) : 0));
    fitted_agd_arm[i] = nprime[i] * pprime[i];

	  for (j in 1:n_int_thin) {
      theta2_bar_cum[(i-1)*n_int_thin + j] = (dot_self(theta_agd_arm_ii[(1 + (i-1)*nint):((i-1)*nint + j*int_thin)]) / (j*int_thin));
    }
  }

}
