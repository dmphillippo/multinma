functions {
#include /include/prior_select.stan
#include /include/count_nonzero.stan

  // -- log-Survival and log-Hazard --
  real lS (row_vector ibasis, real eta, vector scoef) {
    // ibasis = integrated basis evaluated at event/censoring time
    // scoef = row_vector of spline coefficients
    // eta = log rate (linear predictor)
    return - ibasis * scoef .* exp(eta);
  }

  real lh (row_vector basis, real eta, vector scoef) {
    // basis = basis evaluated at event/censoring time
    // scoef = row_vector of spline coefficients
    // eta = log rate (linear predictor)
    return log(basis * scoef) + eta;
  }

  // AgD versions with integration points
  vector lS_a (row_vector ibasis, vector eta, vector scoef) {
    // ibasis = integrated basis evaluated at event/censoring time
    // scoef = row_vector of spline coefficients
    // eta = log rate (linear predictor)
    return - ibasis * scoef * exp(eta);
  }

  vector lh_a (row_vector basis, vector eta, vector scoef) {
    // basis = basis evaluated at event/censoring time
    // scoef = row_vector of spline coefficients
    // eta = log rate (linear predictor)
    return log(basis * scoef) + eta;
  }

  // -- Log likelihood with censoring and truncation --
  real loglik(row_vector time,         // Basis evaluated at event/cens time
              row_vector itime,        // Integrated basis evaluated at event/cens time
              row_vector start_itime,  // Integrated basis evaluated at left interval time
              row_vector delay_itime,  // Integrated basis evaluated at delayed entry time
              int delayed,             // Delayed entry flag (1=delay)
              int status,              // Censoring status
              real eta,                // Linear predictor
              vector scoef) {          // Spline coefficients
    real l;

    if (status == 0) { // Right censored
      l = lS(itime, eta, scoef);
    } else if (status == 1) { // Observed
      l = lS(itime, eta, scoef) + lh(time, eta, scoef);
    } else if (status == 2) { // Left censored
      l = log1m_exp(lS(itime, eta, scoef));
    } else if (status == 3) { // Interval censored
      l = log_diff_exp(lS(start_itime, eta, scoef), lS(itime, eta, scoef));
    }

    // Left truncation
    if (delayed) {
      l -= lS(delay_itime, eta, scoef);
    }

    return l;
  }

  // AgD version with integration points
  vector loglik_a(row_vector time,         // Basis evaluated at event/cens time
                  row_vector itime,        // Integrated basis evaluated at event/cens time
                  row_vector start_itime,  // Integrated basis evaluated at left interval time
                  row_vector delay_itime,  // Integrated basis evaluated at delayed entry time
                  int delayed,             // Delayed entry flag (1=delay)
                  int status,              // Censoring status
                  vector eta,              // Linear predictor
                  vector scoef) {          // Spline coefficients
    vector[num_elements(eta)] l;

    if (status == 0) { // Right censored
      l = lS_a(itime, eta, scoef);
    } else if (status == 1) { // Observed
      l = lS_a(itime, eta, scoef) + lh_a(time, eta, scoef);
    } else if (status == 2) { // Left censored
      l = log1m_exp(lS_a(itime, eta, scoef));
    } else if (status == 3) { // Interval censored
      // l = log_diff_exp(lS_a(start_itime, eta, scoef), lS_a(itime, eta, scoef));
      l = log(exp(lS_a(start_itime, eta, scoef)) - exp(lS_a(itime, eta, scoef)));
    }

    // Left truncation
    if (delayed) {
      l -= lS_a(delay_itime, eta, scoef);
    }

    return l;
  }
}
data {
#include /include/data_common.stan

  // Number of spline coefficients
  int<lower=2> n_scoef;

  // Dirichlet prior on spline coefficients
  real<lower=0> prior_aux_shape;

  // AgD arm-based
  // int<lower=0> narm_agd_arm; // Number of arm-based AgD arms
  int<lower=1> agd_arm_arm[ni_agd_arm]; // Arm indicator for AgD (arm-based) (i.e. picking element of which_RE)

  // Outcomes
  row_vector[n_scoef] ipd_time[ni_ipd];
  row_vector[n_scoef] ipd_itime[ni_ipd];
  row_vector[n_scoef] ipd_start_itime[ni_ipd];
  row_vector[n_scoef] ipd_delay_itime[ni_ipd];
  int<lower=0, upper=1> ipd_delayed[ni_ipd];
  int<lower=0, upper=3> ipd_status[ni_ipd];

  row_vector[n_scoef] agd_arm_time[ni_agd_arm];
  row_vector[n_scoef] agd_arm_itime[ni_agd_arm];
  row_vector[n_scoef] agd_arm_start_itime[ni_agd_arm];
  row_vector[n_scoef] agd_arm_delay_itime[ni_agd_arm];
  int<lower=0, upper=1> agd_arm_delayed[ni_agd_arm];
  int<lower=0, upper=3> agd_arm_status[ni_agd_arm];

  // Study IDs for independent spline coefficients parameters
  int<lower=1> study[narm_ipd + narm_agd_arm];
}
transformed data {
  // Dirichlet prior vector
  vector[n_scoef] prior_aux_shapes = rep_vector(prior_aux_shape, n_scoef);

#include /include/transformed_data_common.stan
}
parameters {
#include /include/parameters_common.stan

  // Spline coefficients
  simplex[n_scoef] scoef[ns_ipd + ns_agd_arm];
}
transformed parameters {
  // Log likelihood contributions
  vector[ni_ipd] log_L_ipd;
  vector[ni_agd_arm] log_L_agd_arm;

#include /include/transformed_parameters_common.stan


  // Evaluate log likelihood
  for (i in 1:ni_ipd) {
    log_L_ipd[i] = loglik(ipd_time[i],
                          ipd_itime[i],
                          ipd_start_itime[i],
                          ipd_delay_itime[i],
                          ipd_delayed[i],
                          ipd_status[i],
                          eta_ipd[i],
                          scoef[study[ipd_arm[i]]]);
  }

  // -- AgD model (arm-based) --
  if (ni_agd_arm) {
    vector[nint_max * ni_agd_arm] eta_agd_arm_noRE = has_offset ?
              X_agd_arm * beta_tilde + offset_agd_arm :
              X_agd_arm * beta_tilde;

    if (nint > 1) { // -- If integration points are used --

      // Integration points are local variable only for efficiency
      vector[nint] eta_agd_arm_ii;
      vector[nint] log_L_ii;

      for (i in 1:ni_agd_arm) {
        eta_agd_arm_ii = eta_agd_arm_noRE[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)];

        // Random effects
        if (RE && which_RE[narm_ipd + agd_arm_arm[i]]) eta_agd_arm_ii += f_delta[which_RE[narm_ipd + agd_arm_arm[i]]];

        // Average likelihood over integration points
        log_L_ii = loglik_a(agd_arm_time[i],
                            agd_arm_itime[i],
                            agd_arm_start_itime[i],
                            agd_arm_delay_itime[i],
                            agd_arm_delayed[i],
                            agd_arm_status[i],
                            eta_agd_arm_ii,
                            scoef[study[narm_ipd + agd_arm_arm[i]]]);

        log_L_agd_arm[i] = log_sum_exp(log_L_ii) - log(nint);
      }
    } else { // -- If no integration --
      for (i in 1:ni_agd_arm) {
        real eta_agd_arm = eta_agd_arm_noRE[i];

        // Random effects
        if (RE && which_RE[narm_ipd + agd_arm_arm[i]]) eta_agd_arm += f_delta[which_RE[narm_ipd + agd_arm_arm[i]]];

        log_L_agd_arm[i] = loglik(agd_arm_time[i],
                                  agd_arm_itime[i],
                                  agd_arm_start_itime[i],
                                  agd_arm_delay_itime[i],
                                  agd_arm_delayed[i],
                                  agd_arm_status[i],
                                  eta_agd_arm,
                                  scoef[study[narm_ipd + agd_arm_arm[i]]]);
      }
    }
  }
}
model {
#include /include/model_common.stan

  // -- Prior on spline coefficients --
  for (s in 1:(ns_ipd + ns_agd_arm)) {
    scoef[s] ~ dirichlet(prior_aux_shapes);
  }

  // -- IPD likelihood --
  target += log_L_ipd;

  // -- AgD likelihood (arm-based) --
  target += log_L_agd_arm;

}
generated quantities {
#include /include/generated_quantities_common.stan

  // Log likelihood
  if (ni_ipd) log_lik[1:ni_ipd] = log_L_ipd;
  if (ni_agd_arm) log_lik[(ni_ipd + 1):(ni_ipd + ni_agd_arm)] = log_L_agd_arm;

  // Residual deviance
  // NOTE: This is actually just the deviance, which is fine for relative fit model comparison
  if (ni_ipd + ni_agd_arm) resdev[1:(ni_ipd + ni_agd_arm)] = -2 * log_lik[1:(ni_ipd + ni_agd_arm)];

  // Fitted values not implemented

}
