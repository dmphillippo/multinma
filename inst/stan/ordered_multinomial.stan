functions {
#include /include/prior_select.stan
#include /include/count_nonzero.stan
}
data {
#include /include/data_common.stan

  // Outcomes
  int<lower=2> ncat;
  int<lower=1, upper=ncat> ipd_r[ni_ipd];
  int<lower=0> agd_arm_r[ni_agd_arm, ncat];

  // Prior on differences between cutpoints
  int<lower=0,upper=5> prior_aux_dist;
  real prior_aux_location;
  real<lower=0> prior_aux_scale;
  real<lower=0> prior_aux_df;

}
transformed data {
  vector[ni_agd_arm] agd_arm_n; // AgD arm sample sizes

#include /include/transformed_data_common.stan

  for (i in 1:ni_agd_arm) {
    agd_arm_n[i] = sum(agd_arm_r[i, ]);
  }
}
parameters {
#include /include/parameters_common.stan

  // Ordered cutoffs on underlying probit-PASI scale
  // "Fixed effect" cutoffs, the same across trials
  positive_ordered[ncat - 2] f_cc;
}
transformed parameters {
  vector[ncat - 1] cc;

  vector[ncat] theta_ipd[ni_ipd]; // IPD transformed predictor

  matrix[nint > 1 ? nint * ni_agd_arm : 0, nint > 1 ? ncat - 1 : 0] theta_agd_arm_ii; // Use these as q_ii intermediates
  vector[ncat] q_agd_arm_bar[ni_agd_arm]; // AgD arm transformed predictor
  vector[ncat] theta_agd_arm_bar[ni_agd_arm]; // AgD arm transformed predictor

#include /include/transformed_parameters_common.stan

  cc[1] = 0;
  cc[2:] = f_cc;

  // -- IPD model --
  // Is this only necessary if link > 2? Since ordered_(logistic|probit) are available
  for (i in 1:ni_ipd) {
    vector[ncat - 1] q_temp;

    // Category 1
    if (link == 1) // logit link
      q_temp[1] = inv_logit(eta_ipd[i] - cc[1]);
    else if (link == 2) // probit link
      q_temp[1] = Phi(eta_ipd[i] - cc[1]);
    else if (link == 3) // cloglog link
      q_temp[1] = inv_cloglog(eta_ipd[i] - cc[1]);

    theta_ipd[i, 1] = 1 - q_temp[1];

    // Categories 2:(ncat - 1)
    for (k in 2:(ncat - 1)) {
      if (link == 1) // logit link
        q_temp[k] = inv_logit(eta_ipd[i] - cc[k]);
      else if (link == 2) // probit link
        q_temp[k] = Phi(eta_ipd[i] - cc[k]);
      else if (link == 3) // cloglog link
        q_temp[k] = inv_cloglog(eta_ipd[i] - cc[k]);

      theta_ipd[i, k] = q_temp[k - 1] - q_temp[k];
    }

    // Category ncat
    theta_ipd[i, ncat] = q_temp[ncat - 1];
  }

  // -- AgD model (arm-based) --
  {
    vector[nint * ni_agd_arm] eta_agd_arm_noRE = has_offset ?
            X_agd_arm * beta_tilde + offset_agd_arm :
            X_agd_arm * beta_tilde;

    if (ni_agd_arm) {
      if (nint > 1) { // -- If integration points are used --

        if (RE) {

          if (link == 1) { // logit link
            for (k in 1:(ncat - 1)) {
              for (i in 1:ni_agd_arm) {
                if (which_RE[narm_ipd + i])
                  theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint), k] = inv_logit(eta_agd_arm_noRE[(1 + (i-1)*nint):(i*nint)] + f_delta[which_RE[narm_ipd + i]] - cc[k]);
                else
                  theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint), k] = inv_logit(eta_agd_arm_noRE[(1 + (i-1)*nint):(i*nint)] - cc[k]);
              }
            }
          } else if (link == 2) { // probit link
            for (k in 1:(ncat - 1)) {
              for (i in 1:ni_agd_arm) {
                if (which_RE[narm_ipd + i])
                  theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint), k] = Phi(eta_agd_arm_noRE[(1 + (i-1)*nint):(i*nint)] + f_delta[which_RE[narm_ipd + i]] - cc[k]);
                else
                  theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint), k] = Phi(eta_agd_arm_noRE[(1 + (i-1)*nint):(i*nint)] - cc[k]);
              }
            }
          } else if (link == 3) { // cloglog link
            for (k in 1:(ncat - 1)) {
              for (i in 1:ni_agd_arm) {
                if (which_RE[narm_ipd + i])
                  theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint), k] = inv_cloglog(eta_agd_arm_noRE[(1 + (i-1)*nint):(i*nint)] + f_delta[which_RE[narm_ipd + i]] - cc[k]);
                else
                  theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint), k] = inv_cloglog(eta_agd_arm_noRE[(1 + (i-1)*nint):(i*nint)] - cc[k]);
              }
            }
          }

        } else {

          if (link == 1) { // logit link
            for (k in 1:(ncat - 1)) {
              for (i in 1:ni_agd_arm) {
                theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint), k] = inv_logit(eta_agd_arm_noRE[(1 + (i-1)*nint):(i*nint)] - cc[k]);
              }
            }
          } else if (link == 2) { // probit link
            for (k in 1:(ncat - 1)) {
              for (i in 1:ni_agd_arm) {
                theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint), k] = Phi(eta_agd_arm_noRE[(1 + (i-1)*nint):(i*nint)] - cc[k]);
              }
            }
          } else if (link == 3) { // cloglog link
            for (k in 1:(ncat - 1)) {
              for (i in 1:ni_agd_arm) {
                theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint), k] = inv_cloglog(eta_agd_arm_noRE[(1 + (i-1)*nint):(i*nint)] - cc[k]);
              }
            }
          }

        }

        for (k in 1:(ncat - 1)) {
          for (i in 1:ni_agd_arm) {
            q_agd_arm_bar[i, k] = mean(theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint), k]);
          }
        }

      } else { // -- If no integration --

        if (RE) {

          if (link == 1) { // logit link
            for (i in 1:ni_agd_arm) {
              for (k in 1:(ncat - 1)) {
                if (which_RE[narm_ipd + i])
                  q_agd_arm_bar[i, k] = inv_logit(eta_agd_arm_noRE[i] + f_delta[which_RE[narm_ipd + i]] - cc[k]);
                else
                  q_agd_arm_bar[i, k] = inv_logit(eta_agd_arm_noRE[i] - cc[k]);
              }
            }
          } else if (link == 2) { // probit link
            for (i in 1:ni_agd_arm) {
              for (k in 1:(ncat - 1)) {
                if (which_RE[narm_ipd + i])
                  q_agd_arm_bar[i, k] = Phi(eta_agd_arm_noRE[i] + f_delta[which_RE[narm_ipd + i]] - cc[k]);
                else
                  q_agd_arm_bar[i, k] = Phi(eta_agd_arm_noRE[i] - cc[k]);
              }
            }
          } else if (link == 3) { // cloglog link
            for (i in 1:ni_agd_arm) {
              for (k in 1:(ncat - 1)) {
                if (which_RE[narm_ipd + i])
                  q_agd_arm_bar[i, k] = inv_cloglog(eta_agd_arm_noRE[i] + f_delta[which_RE[narm_ipd + i]] - cc[k]);
                else
                  q_agd_arm_bar[i, k] = inv_cloglog(eta_agd_arm_noRE[i] - cc[k]);
              }
            }
          }

        } else {

          if (link == 1) { // logit link
            for (i in 1:ni_agd_arm) {
              for (k in 1:(ncat - 1)) {
                q_agd_arm_bar[i, k] = inv_logit(eta_agd_arm_noRE[i] - cc[k]);
              }
            }
          } else if (link == 2) { // probit link
            for (i in 1:ni_agd_arm) {
              for (k in 1:(ncat - 1)) {
                q_agd_arm_bar[i, k] = Phi(eta_agd_arm_noRE[i] - cc[k]);
              }
            }
          } else if (link == 3) { // cloglog link
            for (i in 1:ni_agd_arm) {
              for (k in 1:(ncat - 1)) {
                q_agd_arm_bar[i, k] = inv_cloglog(eta_agd_arm_noRE[i] - cc[k]);
              }
            }
          }

        }
      }


      for (i in 1:ni_agd_arm) {
        // Category 1
        theta_agd_arm_bar[i, 1] = 1 - q_agd_arm_bar[i, 1];

        // Categories 2:(ncat - 1)
        for (k in 2:(ncat - 1))
          theta_agd_arm_bar[i, k] = q_agd_arm_bar[i, k - 1] - q_agd_arm_bar[i, k];

        // Category ncat
        theta_agd_arm_bar[i, ncat] = q_agd_arm_bar[i, ncat - 1];
      }

    }
  }
}
model {
#include /include/model_common.stan

  // -- IPD likelihood --
  for (i in 1:ni_ipd)
    ipd_r[i] ~ categorical(theta_ipd[i]);
  // Use ordered_(logistic|probit) for those links?

  // -- AgD likelihood (arm-based) --
  for (i in 1:ni_agd_arm)
    agd_arm_r[i] ~ multinomial(theta_agd_arm_bar[i]);

  // -- Priors on cutpoints --
  // Implied improper uniform prior on cutpoints if prior_aux_dist = 0
  //   cc ~ uniform(-inf, inf)
  // Otherwise put priors on the differences between cutpoints. Stan will
  // automatically impose the ordering constraints
  if (prior_aux_dist > 0) {
    vector[ncat - 2] diff_cc;
    for (k in 1:(ncat - 2))
      diff_cc[k] = cc[k + 1] - cc[k];
    prior_select_lp(diff_cc, prior_aux_dist, prior_aux_location, prior_aux_scale, prior_aux_df);
  }
}
generated quantities {
  vector[ncat] fitted_ipd[ni_ipd];
  vector[ncat] fitted_agd_arm[ni_agd_arm];
  matrix[ni_agd_arm * n_int_thin, ncat] theta_bar_cum_agd_arm;

#include /include/generated_quantities_common.stan

  // IPD log likelihood and residual deviance
  for (i in 1:ni_ipd) {
    log_lik[i] = categorical_lpmf(ipd_r[i] | theta_ipd[i]);
    resdev[i] = -2 * log_lik[i];
    fitted_ipd[i] = theta_ipd[i];
  }

  // AgD (arm-based) log likelihood and residual deviance
  {
    vector[ncat] dv;
    for (i in 1:ni_agd_arm) {
      log_lik[ni_ipd + i] = multinomial_lpmf(agd_arm_r[i] | theta_agd_arm_bar[i]);
      fitted_agd_arm[i] = agd_arm_n[i] * theta_agd_arm_bar[i];

      // Multinomial residual deviance
      for (k in 1:ncat) {
        dv[k] = agd_arm_r[i, k] == 0 ? 0 : lmultiply(agd_arm_r[i, k], agd_arm_r[i, k] / fitted_agd_arm[i, k]);
      }
      resdev[ni_ipd + i] = 2 * sum(dv);
    }
  }

  // Cumulative integration - note this is for the intermediate q
  for (k in 1:(ncat - 1)) {
    for (i in 1:ni_agd_arm) {
      for (j in 1:n_int_thin) {
        theta_bar_cum_agd_arm[(i - 1)*n_int_thin + j, k] = mean(theta_agd_arm_ii[(1 + (i - 1)*nint):((i - 1)*nint + j*int_thin), k]);
      }
    }
  }

}
