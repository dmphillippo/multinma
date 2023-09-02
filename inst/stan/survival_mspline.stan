functions {
#include /include/prior_select.stan
#include /include/count_nonzero.stan
#include /include/vector_functions.stan

  // -- log-Survival and log-Hazard --
  vector lS (matrix ibasis, vector eta, matrix scoef) {
    // ibasis = integrated basis evaluated at event/censoring time
    // scoef = row_vector of spline coefficients
    // eta = log rate (linear predictor)
    return - rows_dot_product(ibasis, scoef) .* exp(eta);
  }

  vector lh (matrix basis, vector eta, matrix scoef) {
    // basis = basis evaluated at event/censoring time
    // scoef = row_vector of spline coefficients
    // eta = log rate (linear predictor)
    return log(rows_dot_product(basis, scoef)) + eta;
  }

  // AgD versions with integration points
  vector lS_a (row_vector ibasis, vector eta, matrix scoef) {
    // ibasis = integrated basis evaluated at event/censoring time
    // scoef = row_vector of spline coefficients
    // eta = log rate (linear predictor)
    return - scoef * ibasis' .* exp(eta);
  }

  vector lh_a (row_vector basis, vector eta, matrix scoef) {
    // basis = basis evaluated at event/censoring time
    // scoef = row_vector of spline coefficients
    // eta = log rate (linear predictor)
    return log(scoef * basis') + eta;
  }

  // -- Log likelihood with censoring and truncation --
  vector loglik(matrix time,         // Basis evaluated at event/cens time
              matrix itime,        // Integrated basis evaluated at event/cens time
              matrix start_itime,  // Integrated basis evaluated at left interval time
              matrix delay_itime,  // Integrated basis evaluated at delayed entry time
              int[] delayed,             // Delayed entry flag (1=delay)
              int[] status,              // Censoring status
              vector eta,                // Linear predictor
              matrix scoef) {          // Spline coefficients

    vector[num_elements(eta)] l;

    // Right censored
    l = lS(itime, eta, scoef);

    // Observed
    l[which(status, 1)] += lh(time[which(status, 1)], eta[which(status, 1)], scoef);

    // Left censored
    l[which(status, 2)] = log1m_exp(l[which(status, 2)]);

    // Interval censored
    // l = log_diff_exp(lS(start_itime, eta, scoef), lS(itime, eta, scoef));
    l[which(status, 3)] = log(exp(lS(start_itime[which(status, 3)], eta[which(status, 3)], scoef)) - exp(l[which(status, 3)]));

    // Left truncation
    l[which(delayed, 1)] -= lS(delay_itime[which(delayed, 1)], eta[which(delayed, 1)], scoef);

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
                  matrix scoef) {          // Spline coefficients
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

  // Hierarchical logistic prior on spline coefficients
  real<lower=0> prior_aux_a;  // hyperprior on logistic sd sigma ~ Gamma(a, b)
  real<lower=0> prior_aux_b;
  vector[n_scoef-1] prior_aux_mu;  // prior logistic mean

  // AgD arm-based
  // int<lower=0> narm_agd_arm; // Number of arm-based AgD arms
  array[ni_agd_arm] int<lower=1> agd_arm_arm; // Arm indicator for AgD (arm-based) (i.e. picking element of which_RE)

  // Outcomes
  matrix[ni_ipd, n_scoef] ipd_time;
  matrix[ni_ipd, n_scoef] ipd_itime;
  matrix[ni_ipd, n_scoef] ipd_start_itime;
  matrix[ni_ipd, n_scoef] ipd_delay_itime;
  array[ni_ipd] int<lower=0, upper=1> ipd_delayed;
  array[ni_ipd] int<lower=0, upper=3> ipd_status;

  matrix[ni_agd_arm, n_scoef] agd_arm_time;
  matrix[ni_agd_arm, n_scoef] agd_arm_itime;
  matrix[ni_agd_arm, n_scoef] agd_arm_start_itime;
  matrix[ni_agd_arm, n_scoef] agd_arm_delay_itime;
  array[ni_agd_arm] int<lower=0, upper=1> agd_arm_delayed;
  array[ni_agd_arm] int<lower=0, upper=3> agd_arm_status;

  // aux IDs for independent spline coefficients parameters
  // int<lower=0, upper=1> aux_by; // Flag subgroup aux parameters within each arm (1 = yes)
  // int<lower=0, upper=1> aux_reg; // Flag regression on aux pars (1 = yes)
  int<lower=1> aux_id[ni_ipd + ni_agd_arm*nint_max]; //aux_id[ni_ipd + ni_agd_arm*(aux_by || aux_reg ? nint_max : 1)];

  // auxiliary design matrix
  int<lower=0> nX_aux;
  matrix[ni_ipd + nint_max * (ni_agd_arm + ni_agd_contrast), nX_aux] X_aux; // X_aux is Q from QR decomposition if QR = 1
  matrix[QR ? nX_aux : 0, QR ? nX_aux : 0] R_inv_aux;
}
transformed data {
  // Dirichlet prior vector
  // vector[n_scoef] prior_aux_shapes = rep_vector(prior_aux_shape, n_scoef);

  // Number of aux coefficient vectors
  int n_aux = max(aux_id);

  // Split spline Q matrix or X matrix into IPD and AgD rows
  matrix[0, nX_aux] Xauxdummy;
  matrix[ni_ipd, nX_aux] X_aux_ipd = ni_ipd ? X_aux[1:ni_ipd] : Xauxdummy;
  matrix[nint_max * ni_agd_arm, nX_aux] X_aux_agd_arm = ni_agd_arm ? X_aux[(ni_ipd + 1):(ni_ipd + nint_max * ni_agd_arm)] : Xauxdummy;

#include /include/transformed_data_common.stan
}
parameters {
#include /include/parameters_common.stan

  // Spline coefficients
  // simplex[n_scoef] scoef[n_aux];

  // Spline coefficient regression
  matrix[nX_aux, n_scoef-1] beta_aux_tilde;

  // Spline shrinkage SDs and non-centered parameterisation smoothing
  vector<lower=0>[n_aux] sigma;
  vector[n_scoef-1] u_aux[n_aux];
}
transformed parameters {
  // Log likelihood contributions
  vector[ni_ipd] log_L_ipd;
  vector[ni_agd_arm] log_L_agd_arm;

  // Spline coefficients
  matrix[ni_ipd, n_scoef] scoef_ipd;
  matrix[nint_max * ni_agd_arm, n_scoef] scoef_agd_arm;

#include /include/transformed_parameters_common.stan

  // Construct spline coefficients
  if (ni_ipd) {
    matrix[ni_ipd, n_scoef-1] Xb_aux = X_aux_ipd * beta_aux_tilde;
    for (i in 1:ni_ipd) {
      scoef_ipd[i, ] = to_row_vector(softmax(append_row(0, prior_aux_mu + to_vector(Xb_aux[i, ]) + u_aux[aux_id[i]] * sigma[aux_id[i]])));
    }
  }

  if (ni_agd_arm) {
    matrix[ni_agd_arm * nint_max, n_scoef-1] Xb_aux = X_aux_agd_arm * beta_aux_tilde;
    for (i in 1:(ni_agd_arm * nint_max)) {
      scoef_agd_arm[i, ] = to_row_vector(softmax(append_row(0, prior_aux_mu + to_vector(Xb_aux[i, ]) + u_aux[aux_id[ni_ipd + i]] * sigma[aux_id[ni_ipd + i]])));
    }
  }


  // Evaluate log likelihood
  if (ni_ipd) {
    log_L_ipd = loglik(ipd_time,
                       ipd_itime,
                       ipd_start_itime,
                       ipd_delay_itime,
                       ipd_delayed,
                       ipd_status,
                       eta_ipd,
                       scoef_ipd);
  }

  // -- AgD model (arm-based) --
  if (ni_agd_arm) {
    vector[nint_max * ni_agd_arm] eta_agd_arm_noRE = has_offset ?
              X_agd_arm * beta_tilde + offset_agd_arm :
              X_agd_arm * beta_tilde;

    if (nint_max > 1) { // -- If integration points are used --

      // Integration points are local variable only for efficiency
      vector[nint] eta_agd_arm_ii;
      vector[nint] log_L_ii;

      for (i in 1:ni_agd_arm) {
        eta_agd_arm_ii = eta_agd_arm_noRE[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)];

        // Random effects
        if (RE && which_RE[narm_ipd + agd_arm_arm[i]]) eta_agd_arm_ii += f_delta[which_RE[narm_ipd + agd_arm_arm[i]]];

        // Average likelihood over integration points
        // if (aux_by) for (j in 1:nint) {
          log_L_ii = loglik_a(agd_arm_time[i],
                               agd_arm_itime[i],
                               agd_arm_start_itime[i],
                               agd_arm_delay_itime[i],
                               agd_arm_delayed[i],
                               agd_arm_status[i],
                               eta_agd_arm_ii,
                               scoef_agd_arm[(1 + (i-1)*nint_max):((i-1)*nint_max + nint), ]);
        // } else {
        //   log_L_ii = loglik_a(agd_arm_time[i],
        //                       agd_arm_itime[i],
        //                       agd_arm_start_itime[i],
        //                       agd_arm_delay_itime[i],
        //                       agd_arm_delayed[i],
        //                       agd_arm_status[i],
        //                       eta_agd_arm_ii,
        //                       scoef_agd_arm[aux_id[ni_ipd + i]]);
        // }

        log_L_agd_arm[i] = log_sum_exp(log_L_ii) - log(nint);
      }
    } else { // -- If no integration --
      // for (i in 1:ni_agd_arm) {
        vector[ni_agd_arm] eta_agd_arm = eta_agd_arm_noRE;

        // Random effects
        if (RE) for (i in 1:ni_agd_arm) {
          if (which_RE[narm_ipd + agd_arm_arm[i]]) eta_agd_arm[i] += f_delta[which_RE[narm_ipd + agd_arm_arm[i]]];
        }

        log_L_agd_arm = loglik(agd_arm_time,
                                  agd_arm_itime,
                                  agd_arm_start_itime,
                                  agd_arm_delay_itime,
                                  agd_arm_delayed,
                                  agd_arm_status,
                                  eta_agd_arm,
                                  scoef_agd_arm);
      // }
    }
  }
}
model {
#include /include/model_common.stan

  // -- Prior on spline coefficients --
  for (i in 1:n_aux) u_aux[i] ~ logistic(0, 1);
  sigma ~ gamma(prior_aux_a, prior_aux_b);
  for (i in 1:(n_scoef-1)) prior_select_lp(beta_aux_tilde[, i], prior_reg_dist, prior_reg_location, prior_reg_scale, prior_reg_df);

  // -- IPD likelihood --
  target += log_L_ipd;

  // -- AgD likelihood (arm-based) --
  target += log_L_agd_arm;

}
generated quantities {
  // Baseline spline coefficients
  vector[n_scoef] scoef[n_aux];

  // Regression on spline coefficients
  matrix[nX_aux, n_scoef-1] beta_aux = QR ? R_inv_aux * beta_aux_tilde : beta_aux_tilde;

#include /include/generated_quantities_common.stan

  // Baseline spline coefficients
  for (i in 1:n_aux) {
    scoef[i] = softmax(append_row(0, prior_aux_mu + u_aux[i] * sigma[i]));
  }

  // Log likelihood
  if (ni_ipd) log_lik[1:ni_ipd] = log_L_ipd;
  if (ni_agd_arm) log_lik[(ni_ipd + 1):(ni_ipd + ni_agd_arm)] = log_L_agd_arm;

  // Residual deviance
  // NOTE: This is actually just the deviance, which is fine for relative fit model comparison
  if (ni_ipd + ni_agd_arm) resdev[1:(ni_ipd + ni_agd_arm)] = -2 * log_lik[1:(ni_ipd + ni_agd_arm)];

  // Fitted values not implemented

}
