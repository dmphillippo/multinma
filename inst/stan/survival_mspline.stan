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

  // Versions with a single scoef vector for all individuals
  vector lS2 (matrix ibasis, vector eta, vector scoef) {
    // ibasis = integrated basis evaluated at event/censoring time
    // scoef = row_vector of spline coefficients
    // eta = log rate (linear predictor)
    return - (ibasis * scoef) .* exp(eta);
  }

  vector lh2 (matrix basis, vector eta, vector scoef) {
    // basis = basis evaluated at event/censoring time
    // scoef = row_vector of spline coefficients
    // eta = log rate (linear predictor)
    return log(basis * scoef) + eta;
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

  matrix lS_a2 (matrix ibasis, matrix eta, vector scoef) {
    // ibasis = integrated basis evaluated at event/censoring time
    // scoef = row_vector of spline coefficients
    // eta = log rate (linear predictor)
    int n = rows(ibasis);
    vector[n] ibs = - ibasis * scoef;
    matrix[rows(eta), n] eeta = exp(eta);
    matrix[rows(eta), n] l;
    for (i in 1:n) l[,i] = ibs[i] * eeta[,i];
    return l;
  }

  matrix lh_a2 (matrix basis, matrix eta, vector scoef) {
    // basis = basis evaluated at event/censoring time
    // scoef = row_vector of spline coefficients
    // eta = log rate (linear predictor)
    int n = rows(basis);
    vector[n] lbs = log(basis * scoef);
    matrix[rows(eta), n] l;
    for (i in 1:n) l[,i] = lbs[i] + eta[,i];
    return l;
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
    array[3] int nw = nwhich_all(status, 3);
    int nd = nwhich(delayed, 1);

    // Right censored
    l = lS(itime, eta, scoef);

    // Observed
    if (nw[1]) {
      array[nw[1]] int w1 = which(status, 1);
      l[w1] += lh(time[w1], eta[w1], scoef[w1]);
    }

    // Left censored
    if (nw[2]) {
      array[nw[2]] int w2 = which(status, 2);
      l[w2] = log1m_exp(l[w2]);
    }

    // Interval censored
    if (nw[3]) {
      array[nw[3]] int w3 = which(status, 3);
      // l[w3] = log(exp(lS(start_itime[w3], eta[w3], scoef[w3])) - exp(l[w3]));
      l[w3] = log_diff_exp(lS(start_itime[w3], eta[w3], scoef[w3]), l[w3]);
    }

    // Left truncation
    if (nd) {
      array[nd] int wd = which(delayed, 1);
      l[wd] -= lS(delay_itime[wd], eta[wd], scoef[wd]);
    }

    return l;
  }

  // Version with single scoef vector for all individuals
  vector loglik2(matrix time,         // Basis evaluated at event/cens time
              matrix itime,        // Integrated basis evaluated at event/cens time
              matrix start_itime,  // Integrated basis evaluated at left interval time
              matrix delay_itime,  // Integrated basis evaluated at delayed entry time
              int[] delayed,             // Delayed entry flag (1=delay)
              int[] status,              // Censoring status
              vector eta,                // Linear predictor
              vector scoef) {          // Spline coefficients

    vector[num_elements(eta)] l;
    array[3] int nw = nwhich_all(status, 3);
    int nd = nwhich(delayed, 1);

    // Right censored
    l = lS2(itime, eta, scoef);

    // Observed
    if (nw[1]) {
      array[nw[1]] int w1 = which(status, 1);
      l[w1] += lh2(time[w1], eta[w1], scoef);
    }

    // Left censored
    if (nw[2]) {
      array[nw[2]] int w2 = which(status, 2);
      l[w2] = log1m_exp(l[w2]);
    }

    // Interval censored
    if (nw[3]) {
      array[nw[3]] int w3 = which(status, 3);
      // l[w3] = log(exp(lS2(start_itime[w3], eta[w3], scoef)) - exp(l[w3]));
      l[w3] = log_diff_exp(lS2(start_itime[w3], eta[w3], scoef), l[w3]);
    }

    // Left truncation
    if (nd) {
      array[nd] int wd = which(delayed, 1);
      l[wd] -= lS2(delay_itime[wd], eta[wd], scoef);
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
                  matrix scoef) {          // Spline coefficients
    vector[num_elements(eta)] l;

    if (status == 0) { // Right censored
      l = lS_a(itime, eta, scoef);
    } else if (status == 1) { // Observed
      l = lS_a(itime, eta, scoef) + lh_a(time, eta, scoef);
    } else if (status == 2) { // Left censored
      l = log1m_exp(lS_a(itime, eta, scoef));
    } else if (status == 3) { // Interval censored
      l = log_diff_exp(lS_a(start_itime, eta, scoef), lS_a(itime, eta, scoef));
      // l = log(exp(lS_a(start_itime, eta, scoef)) - exp(lS_a(itime, eta, scoef)));
    }

    // Left truncation
    if (delayed) {
      l -= lS_a(delay_itime, eta, scoef);
    }

    return l;
  }

  // AgD version with single scoef vector for all individuals
  matrix loglik_a2(matrix time,         // Basis evaluated at event/cens time
              matrix itime,        // Integrated basis evaluated at event/cens time
              matrix start_itime,  // Integrated basis evaluated at left interval time
              matrix delay_itime,  // Integrated basis evaluated at delayed entry time
              int[] delayed,             // Delayed entry flag (1=delay)
              int[] status,              // Censoring status
              matrix eta,                // Linear predictor
              vector scoef) {          // Spline coefficients

    matrix[rows(eta), cols(eta)] l;
    array[3] int nw = nwhich_all(status, 3);
    int nd = nwhich(delayed, 1);

    // Right censored
    l = lS_a2(itime, eta, scoef);

    // Observed
    if (nw[1]) {
      array[nw[1]] int w1 = which(status, 1);
      l[,w1] += lh_a2(time[w1], eta[,w1], scoef);
    }

    // Left censored
    if (nw[2]) {
      array[nw[2]] int w2 = which(status, 2);
      l[,w2] = log1m_exp(l[,w2]);
    }

    // Interval censored
    if (nw[3]) {
      array[nw[3]] int w3 = which(status, 3);
      // l[,w3] = log(exp(lS_a2(start_itime[w3], eta[,w3], scoef)) - exp(l[,w3]));
      l[,w3] = log_diff_exp(lS_a2(start_itime[w3], eta[,w3], scoef), l[,w3]);
    }

    // Left truncation
    if (nd) {
      array[nd] int wd = which(delayed, 1);
      l[,wd] -= lS_a2(delay_itime[wd], eta[,wd], scoef);
    }

    return l;
  }

}
data {
#include /include/data_common.stan

  // Number of spline coefficients
  int<lower=2> n_scoef;

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
  int<lower=0, upper=1> aux_int; // Flag aux regression needs integration (1 = yes)
  array[ni_ipd + ni_agd_arm*(aux_int ? nint_max : 1)] int<lower=1> aux_id;

  // auxiliary design matrix
  int<lower=0> nX_aux;
  matrix[ni_ipd + (aux_int ? nint_max : 1) * (ni_agd_arm + ni_agd_contrast), nX_aux] X_aux; // X_aux is Q from QR decomposition if QR = 1
  matrix[QR ? nX_aux : 0, QR ? nX_aux : 0] R_inv_aux;

  // Hierarchical logistic prior on spline coefficients
  int<lower=0,upper=6> prior_hyper_dist;
  real<lower=0> prior_hyper_location;
  real<lower=0> prior_hyper_scale;
  real<lower=0> prior_hyper_df;
  // real<lower=0> prior_hyper_shape;
  array[nX_aux ? 1 : max(aux_id)] vector[n_scoef-1] prior_aux_location;  // prior logistic mean

  // Hyperprior on aux regression smooths
  int<lower=0,upper=6> prior_reg_hyper_dist;
  real<lower=0> prior_reg_hyper_location;
  real<lower=0> prior_reg_hyper_scale;
  real<lower=0> prior_reg_hyper_df;
}
transformed data {
  // Number of aux coefficient vectors
  int n_aux = max(aux_id);

  // Scaling for spline random walk SD sigma
  real sigma_scale = sqrt(n_scoef-1);

  // Split aux_id by IPD and AgD rows
  array[ni_ipd] int<lower=1> aux_id_ipd = aux_id[1:ni_ipd];
  array[ni_agd_arm*(aux_int ? nint_max : 1)] int<lower=1> aux_id_agd_arm = aux_id[(ni_ipd + 1):(ni_ipd + ni_agd_arm*(aux_int ? nint_max : 1))];

  // Aux ID indexing arrays
  array[aux_int == 0 ? n_aux : 0] int ni_aux_id_ipd = nwhich_all(aux_id_ipd, n_aux);
  array[aux_int == 0 ? n_aux : 0, max(ni_aux_id_ipd)] int wi_aux_id_ipd;
  array[aux_int == 0 ? n_aux : 0] int ni_aux_id_agd_arm = nwhich_all(aux_id_agd_arm, n_aux);
  array[aux_int == 0 ? n_aux : 0, max(ni_aux_id_agd_arm)] int wi_aux_id_agd_arm;

  // Split spline Q matrix or X matrix into IPD and AgD rows
  matrix[0, nX_aux] Xauxdummy;
  matrix[ni_ipd, nX_aux] X_aux_ipd = ni_ipd ? X_aux[1:ni_ipd] : Xauxdummy;
  matrix[(aux_int ? nint_max : 1) * ni_agd_arm, nX_aux] X_aux_agd_arm = ni_agd_arm ? X_aux[(ni_ipd + 1):(ni_ipd + (aux_int ? nint_max : 1) * ni_agd_arm)] : Xauxdummy;

#include /include/transformed_data_common.stan

  if (aux_int == 0) for (i in 1:n_aux) {
    if (ni_aux_id_ipd[i]) wi_aux_id_ipd[i, 1:ni_aux_id_ipd[i]] = which(aux_id_ipd, i);
    if (ni_aux_id_agd_arm[i]) wi_aux_id_agd_arm[i, 1:ni_aux_id_agd_arm[i]] = which(aux_id_agd_arm, i);
  }
}
parameters {
#include /include/parameters_common.stan

  // Spline coefficient regression with shrinkage over time
  matrix[nX_aux, n_scoef-1] u_beta_aux;
  vector<lower=0>[nX_aux] sigma_beta;

  // Spline shrinkage SDs and non-centered parameterisation smoothing
  vector<lower=0>[n_aux] sigma;
  array[n_aux] vector[n_scoef-1] u_aux;
}
transformed parameters {
  // Log likelihood contributions
  vector[ni_ipd] log_L_ipd;
  vector[ni_agd_arm] log_L_agd_arm;

  // Spline coefficients
  array[n_aux] vector[n_scoef-1] lscoef; // logit baseline spline coefficients
  array[nX_aux ? 0 : n_aux] vector[n_scoef] scoef_temp;

  // Spline coefficient regression
  matrix[nX_aux, n_scoef-1] beta_aux;

#include /include/transformed_parameters_common.stan

  // Shrinkage regression on aux pars
  for (i in 1:nX_aux) beta_aux[i, ] = u_beta_aux[i, ] * sigma_beta[i];

  // Construct spline coefficients with random walk prior around constant hazard
  if (nX_aux) {
    for (i in 1:n_aux) {
      lscoef[i] = cumulative_sum(u_aux[i]) * sigma[i] / sigma_scale + prior_aux_location[1];
    }
  } else {
    for (i in 1:n_aux) {
      lscoef[i] = cumulative_sum(u_aux[i]) * sigma[i] / sigma_scale + prior_aux_location[i];
      scoef_temp[i] = softmax(append_row(0, lscoef[i]));
    }
  }


  // Evaluate log likelihood
  if (ni_ipd) {

    if (aux_int) {
      matrix[ni_ipd, n_scoef] scoef_ipd;
      matrix[ni_ipd, n_scoef-1] Xb_aux = X_aux_ipd * beta_aux;
      for (i in 1:ni_ipd) {
        scoef_ipd[i, ] = to_row_vector(softmax(append_row(0, lscoef[aux_id_ipd[i]] + to_vector(Xb_aux[i, ]))));
      }

      log_L_ipd = loglik(ipd_time,
                         ipd_itime,
                         ipd_start_itime,
                         ipd_delay_itime,
                         ipd_delayed,
                         ipd_status,
                         eta_ipd,
                         scoef_ipd);
    } else {
      // Loop over aux parameters (i.e. by study, possibly further stratified by aux_by)
      for (i in 1:n_aux) {
        int ni = ni_aux_id_ipd[i];

        if (ni) {
          array[ni] int wi = wi_aux_id_ipd[i, 1:ni];

          if (nX_aux) {
            row_vector[n_scoef-1] Xb_auxi = X_aux_ipd[first(aux_id_ipd, i), ] * beta_aux;
            log_L_ipd[wi] = loglik2(ipd_time[wi],
                               ipd_itime[wi],
                               ipd_start_itime[wi],
                               ipd_delay_itime[wi],
                               ipd_delayed[wi],
                               ipd_status[wi],
                               eta_ipd[wi],
                               softmax(append_row(0, lscoef[i] + to_vector(Xb_auxi))));
          } else {
            log_L_ipd[wi] = loglik2(ipd_time[wi],
                               ipd_itime[wi],
                               ipd_start_itime[wi],
                               ipd_delay_itime[wi],
                               ipd_delayed[wi],
                               ipd_status[wi],
                               eta_ipd[wi],
                               scoef_temp[i]);
          }
        }
      }
    }
  }

  // -- AgD model (arm-based) --
  if (ni_agd_arm) {
    vector[nint_max * ni_agd_arm] eta_agd_arm_noRE = has_offset ?
              X_agd_arm * beta_tilde + offset_agd_arm :
              X_agd_arm * beta_tilde;


    if (nint_max > 1) { // -- If integration points are used --

      if (aux_int) {
        for (i in 1:ni_agd_arm) {

          vector[nint] eta_agd_arm_ii;
          vector[nint] log_L_ii;
          matrix[nint, n_scoef-1] Xb_aux = X_aux_agd_arm[(1 + (i-1)*nint_max):((i-1)*nint_max + nint), ] * beta_aux;
          matrix[nint, n_scoef] scoef_agd_arm;

          eta_agd_arm_ii = eta_agd_arm_noRE[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)];

          // Random effects
          if (RE && which_RE[narm_ipd + agd_arm_arm[i]]) eta_agd_arm_ii += f_delta[which_RE[narm_ipd + agd_arm_arm[i]]];

          // Average likelihood over integration points
          for (j in 1:nint) {
            scoef_agd_arm[j, ] = to_row_vector(softmax(append_row(0, lscoef[aux_id_agd_arm[(i-1)*nint_max + j]] + to_vector(Xb_aux[j, ]))));
          }
          log_L_ii = loglik_a(agd_arm_time[i],
                               agd_arm_itime[i],
                               agd_arm_start_itime[i],
                               agd_arm_delay_itime[i],
                               agd_arm_delayed[i],
                               agd_arm_status[i],
                               eta_agd_arm_ii,
                               scoef_agd_arm);

          log_L_agd_arm[i] = log_sum_exp(log_L_ii) - log(nint);

        }
      } else {
        // Loop over aux parameters first (i.e. by study, possibly further stratified by aux_by) for efficiency
        for (i in 1:n_aux) {
          int ni = ni_aux_id_agd_arm[i];

          if (ni) {
            array[ni] int wi = wi_aux_id_agd_arm[i, 1:ni];
            matrix[nint, ni] eta_agd_arm_iim;
            matrix[nint, ni] log_L_iim;

            for (j in 1:ni) {
              eta_agd_arm_iim[, j] = eta_agd_arm_noRE[(1 + (wi[j]-1)*nint_max):((wi[j]-1)*nint_max + nint)];
              if (RE && which_RE[narm_ipd + agd_arm_arm[wi[j]]]) eta_agd_arm_iim[, j] += f_delta[which_RE[narm_ipd + agd_arm_arm[wi[j]]]];
            }

            if (nX_aux) {
              row_vector[n_scoef-1] Xb_auxi = X_aux_agd_arm[first(aux_id_agd_arm, i), ] * beta_aux;
              log_L_iim = loglik_a2(agd_arm_time[wi],
                                   agd_arm_itime[wi],
                                   agd_arm_start_itime[wi],
                                   agd_arm_delay_itime[wi],
                                   agd_arm_delayed[wi],
                                   agd_arm_status[wi],
                                   eta_agd_arm_iim,
                                   softmax(append_row(0, lscoef[i] + to_vector(Xb_auxi))));
            } else {
              log_L_iim = loglik_a2(agd_arm_time[wi],
                                  agd_arm_itime[wi],
                                  agd_arm_start_itime[wi],
                                  agd_arm_delay_itime[wi],
                                  agd_arm_delayed[wi],
                                  agd_arm_status[wi],
                                  eta_agd_arm_iim,
                                  scoef_temp[i]);
            }

            for (j in 1:ni) {
              log_L_agd_arm[wi[j]] = log_sum_exp(log_L_iim[, j]) - log(nint);
            }
          }
        }
      }

    } else { // -- If no integration --
      vector[ni_agd_arm] eta_agd_arm = eta_agd_arm_noRE;

      // Random effects
      if (RE) for (i in 1:ni_agd_arm) {
        if (which_RE[narm_ipd + agd_arm_arm[i]]) eta_agd_arm[i] += f_delta[which_RE[narm_ipd + agd_arm_arm[i]]];
      }

      if (aux_int) {
        matrix[ni_agd_arm, n_scoef] scoef_agd_arm;
        matrix[ni_agd_arm, n_scoef-1] Xb_aux = X_aux_agd_arm * beta_aux;
        for (i in 1:ni_agd_arm) {
          scoef_agd_arm[i, ] = to_row_vector(softmax(append_row(0, lscoef[aux_id_agd_arm[i]] + to_vector(Xb_aux[i, ]))));
        }

        log_L_agd_arm = loglik(agd_arm_time,
                               agd_arm_itime,
                               agd_arm_start_itime,
                               agd_arm_delay_itime,
                               agd_arm_delayed,
                               agd_arm_status,
                               eta_agd_arm,
                               scoef_agd_arm);
      } else {
        // Loop over aux parameters (i.e. by study, possibly further stratified by aux_by)
        for (i in 1:n_aux) {
          int ni = ni_aux_id_agd_arm[i];

          if (ni) {
            array[ni] int wi = wi_aux_id_agd_arm[i, 1:ni];

            if (nX_aux) {
              row_vector[n_scoef-1] Xb_auxi = X_aux_agd_arm[first(aux_id_agd_arm, i), ] * beta_aux;
              log_L_agd_arm[wi] = loglik2(agd_arm_time[wi],
                                          agd_arm_itime[wi],
                                          agd_arm_start_itime[wi],
                                          agd_arm_delay_itime[wi],
                                          agd_arm_delayed[wi],
                                          agd_arm_status[wi],
                                          eta_agd_arm[wi],
                                          softmax(append_row(0, lscoef[i] + to_vector(Xb_auxi))));
            } else {
              log_L_agd_arm[wi] = loglik2(agd_arm_time[wi],
                                          agd_arm_itime[wi],
                                          agd_arm_start_itime[wi],
                                          agd_arm_delay_itime[wi],
                                          agd_arm_delayed[wi],
                                          agd_arm_status[wi],
                                          eta_agd_arm[wi],
                                          scoef_temp[i]);
            }
          }
        }
      }
    }
  }
}
model {
#include /include/model_common.stan

  // -- Prior on spline coefficients --
  for (i in 1:n_aux) u_aux[i] ~ logistic(0, 1);

  // -- Hyperprior on spline sd --
    prior_select_lp(sigma, prior_hyper_dist, prior_hyper_location, prior_hyper_scale, prior_hyper_df);

  // -- Smoothing prior on aux regression beta --
  for (i in 1:(n_scoef -1)) u_beta_aux[, i] ~ std_normal();
  prior_select_lp(sigma_beta, prior_reg_hyper_dist, prior_reg_hyper_location, prior_reg_hyper_scale, prior_reg_hyper_df);

  // -- IPD likelihood --
  target += log_L_ipd;

  // -- AgD likelihood (arm-based) --
  target += log_L_agd_arm;

}
generated quantities {
  // Baseline spline coefficients
  array[n_aux] vector[n_scoef] scoef;

  // Regression on spline coefficients
  // matrix[nX_aux, n_scoef-1] beta_aux = QR ? R_inv_aux * beta_aux_tilde : beta_aux_tilde;

#include /include/generated_quantities_common.stan

  // Baseline spline coefficients
  for (i in 1:n_aux) {
    scoef[i] = nX_aux ? softmax(append_row(0, lscoef[i])) : scoef_temp[i];
  }

  // Log likelihood
  if (ni_ipd) log_lik[1:ni_ipd] = log_L_ipd;
  if (ni_agd_arm) log_lik[(ni_ipd + 1):(ni_ipd + ni_agd_arm)] = log_L_agd_arm;

  // Residual deviance
  // NOTE: This is actually just the deviance, which is fine for relative fit model comparison
  if (ni_ipd + ni_agd_arm) resdev[1:(ni_ipd + ni_agd_arm)] = -2 * log_lik[1:(ni_ipd + ni_agd_arm)];

  // Fitted values not implemented

}
