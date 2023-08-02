functions {
#include /include/prior_select.stan
#include /include/count_nonzero.stan
#include /include/vector_functions.stan

  // -- Generalised Gamma ldpf --
  real gengamma_lpdf(real y, real mu, real sigma, real k) {
    // Parameterisation of Lawless
	  real Q = pow(k, -0.5);
	  real w = Q * (log(y) - mu) / sigma;
	  return -log(sigma) - log(y) - 0.5 * log(k) * (1 - 2 * k) + k * (w - exp(w)) - lgamma(k);
  }

  //-- log Survival functions --
  vector lS(int dist, vector y, vector eta, vector aux, vector aux2) {
    // aux is shape, except for lognormal sdlog, gengamma sigma
    // aux2 is only for gengamma k

    int n = num_elements(y);
    vector[n] out;

    if (dist == 1) { // Exponential
      out = -y .* exp(eta);
    } else if (dist == 2) { // Weibull
      out = -pow_vec(y, aux) .* exp(eta);
    } else if (dist == 3) { // Gompertz
      out = -exp(eta)./aux .* expm1(aux .* y);
    } else if (dist == 4) { // Exponential AFT
      out = -y .* exp(-eta);
    } else if (dist == 5) { // Weibull AFT
      out = -pow_vec(y, aux) .* exp(-aux .* eta);
    } else if (dist == 6) { // log Normal
      // aux is sdlog
      // out = log1m(Phi((log(y) - eta) ./ aux));
      for (i in 1:n) out[i] = lognormal_lccdf(y[i] | eta[i], aux[i]);
    } else if (dist == 7) { // log logistic
      out = -log1p(pow_vec(y ./ exp(eta), aux));
    } else if (dist == 8) { // Gamma
      vector[n] eeta = exp(-eta);
      for (i in 1:n) out[i] = gamma_lccdf(y[i] | aux[i], eeta[i]);
    } else if (dist == 9) { // Generalised Gamma
      vector[n] Q = inv(sqrt(aux2));
      vector[n] w = exp(Q .* (log(y) - eta) ./ aux) .* aux2;
      for (i in 1:n) out[i] = log1m(gamma_p(aux2[i], w[i]));
    }

    return out;
  }

  // AgD version with integration points
  vector lS_a(int dist, real y, vector eta, real aux, real aux2) {
    // aux is shape, except for lognormal sdlog, gengamma sigma
    // aux2 is only for gengamma k

    int n = num_elements(eta);
    vector[n] out;

    if (dist == 1) { // Exponential
      out = -y * exp(eta);
    } else if (dist == 2) { // Weibull
      out = -pow(y, aux) * exp(eta);
    } else if (dist == 3) { // Gompertz
      out = -exp(eta)/aux * expm1(aux * y);
    } else if (dist == 4) { // Exponential AFT
      out = -y * exp(-eta);
    } else if (dist == 5) { // Weibull AFT
      out = -pow(y, aux) * exp(-aux * eta);
    } else if (dist == 6) { // log Normal
      // aux is sdlog
      // lS = log1m(Phi((log(y) - eta) / aux));
      for (i in 1:n) out[i] = lognormal_lccdf(y | eta[i], aux);
    } else if (dist == 7) { // log logistic
      out = -log1p(pow_vec2(y ./ exp(eta), aux));
    } else if (dist == 8) { // Gamma
      vector[n] eeta = exp(-eta);
      for (i in 1:n) out[i] = gamma_lccdf(y | aux, eeta[i]);
    } else if (dist == 9) { // Generalised Gamma
      real Q = pow(aux2, -0.5);
      vector[n] w = exp(Q * (log(y) - eta) / aux) * aux2;
      for (i in 1:n) out[i] = log1m(gamma_p(aux2, w[i]));
    }

    return out;
  }

  // -- log Hazard functions --
  vector lh(int dist, vector y, vector eta, vector aux, vector aux2) {
    // aux is shape, except for lognormal sdlog, gengamma scale
    // aux2 is only for gengamma shape

    int n = num_elements(y);
    vector[n] out;

    if (dist == 1) { // Exponential
      out = eta;
    } else if (dist == 2) { // Weibull
      out = log(aux) + eta + lmultiply_vec(aux - 1, y);
    } else if (dist == 3) { // Gompertz
      out = eta + (aux .* y);
    } else if (dist == 4) { // Exponential AFT
      out = -eta;
    } else if (dist == 5) { // Weibull AFT
      out = log(aux) - (aux .* eta) + lmultiply_vec(aux - 1, y);
    } else if (dist == 6) { // log Normal
      // aux is sdlog
      // out = lognormal_lpdf(y | eta, aux) - log1m(Phi((log(y) - eta) ./ aux));
      for (i in 1:n) out[i] = lognormal_lpdf(y[i] | eta[i], aux[i]) - lognormal_lccdf(y[i] | eta[i], aux[i]);
    } else if (dist == 7) { // log logistic
      out = log(aux) - eta + (aux - 1).*(log(y) - eta) - log1p(pow_vec(y ./ exp(eta), aux));
    } else if (dist == 8) { // Gamma
      vector[n] eeta = exp(-eta);
      for (i in 1:n) out[i] = gamma_lpdf(y[i] | aux[i], eeta[i]) - gamma_lccdf(y[i] | aux[i], eeta[i]);
    } else if (dist == 9) { // Generalised Gamma
      // Not used, lpdf used directly
    }

    return out;
  }

  // AgD version with integration points
  vector lh_a(int dist, real y, vector eta, real aux, real aux2) {
    // aux is shape, except for lognormal sdlog, gengamma scale
    // aux2 is only for gengamma shape

    int n = num_elements(eta);
    vector[n] out;

    if (dist == 1) { // Exponential
      out = eta;
    } else if (dist == 2) { // Weibull
      out = log(aux) + eta + lmultiply(aux - 1, y);
    } else if (dist == 3) { // Gompertz
      out = eta + (aux * y);
    } else if (dist == 4) { // Exponential AFT
      out = -eta;
    } else if (dist == 5) { // Weibull AFT
      out = log(aux) - (aux * eta) + lmultiply(aux - 1, y);
    } else if (dist == 6) { // log Normal
      // aux is sdlog
      // out = lognormal_lpdf(y | eta, aux) - log1m(Phi((log(y) - eta) ./ aux));
      for (i in 1:n) out[i] = lognormal_lpdf(y| eta[i], aux) - lognormal_lccdf(y | eta[i], aux);
    } else if (dist == 7) { // log logistic
      out = log(aux) - eta + (aux - 1)*(log(y) - eta) - log1p(pow_vec2(y ./ exp(eta), aux));
    } else if (dist == 8) { // Gamma
      vector[n] eeta = exp(-eta);
      for (i in 1:n) out[i] = gamma_lpdf(y | aux, eeta[i]) - gamma_lccdf(y | aux, eeta[i]);
    } else if (dist == 9) { // Generalised Gamma
      // Not used, lpdf used directly
    }

    return out;
  }

  // -- Log likelihood with censoring and truncation --
  vector loglik(int dist, vector time, vector start_time, vector delay_time, int[] status, vector eta, vector aux, vector aux2) {
    int n = num_elements(eta);
    vector[n] l;

    // Status indicators
    int nwhich0 = (dist == 6 || dist == 8 || dist == 9) ? num_elements(which(status, 0)) : 0;
    int nwhich1 = !(dist == 6 || dist == 8 || dist == 9) ? num_elements(which(status, 1)) : 0;
    int nwhich2 = num_elements(which(status, 2));
    int nwhich3 = num_elements(which(status, 3));
    int nwhichd = num_elements(which_gt0(delay_time));
    int which0[nwhich0];
    int which1[nwhich1];
    int which2[nwhich2];
    int which3[nwhich3];
    int whichd[nwhichd];
    if (nwhich0) which0 = which(status, 0);
    if (nwhich1) which1 = which(status, 1);
    if (nwhich2) which2 = which(status, 2);
    if (nwhich3) which3 = which(status, 3);
    if (nwhichd) whichd = which_gt0(delay_time);

    // Make certain models more efficient by using ldpf directly
    if (dist == 6 || dist == 8 || dist == 9) {
      // Right censored
      l[which0] = lS(dist, time[which0], eta[which0], aux[which0], aux2[which0]);

      // Observed
      if (dist == 6) for (i in 1:n) if (status[i] == 1) l[i] = lognormal_lpdf(time[i] | eta[i], aux[i]);
      if (dist == 8) for (i in 1:n) if (status[i] == 1) l[i] = gamma_lpdf(time[i] | aux[i], exp(-eta[i]));
      if (dist == 9) for (i in 1:n) if (status[i] == 1) l[i] = gengamma_lpdf(time[i] | eta[i], aux[i], aux2[i]);

    } else {

      // Right censored
      l = lS(dist, time, eta, aux, aux2);

      // Observed
      l[which1] += lh(dist, time[which1], eta[which1], aux[which1], aux2[which1]);
    }

    // Left censored
    l[which2] = log1m_exp(l[which2]);

    // Interval censored
    // l = log_diff_exp(lS(dist, start_time, eta, aux, aux2), lS(dist, time, eta, aux, aux2));
    l[which3] = log(exp(lS(dist, start_time[which3], eta[which3], aux[which3], aux2[which3])) - exp(l[which3]));

    // Left truncation
    l[whichd] -= lS(dist, delay_time[whichd], eta[whichd], aux[whichd], aux2[whichd]);

    return l;
  }

  vector loglik_a(int dist, real time, real start_time, real delay_time, int status, vector eta, real aux, real aux2) {
    int n = num_elements(eta);
    vector[n] l;

    if (status == 0) { // Right censored
      l = lS_a(dist, time, eta, aux, aux2);
    } else if (status == 1) { // Observed
      // Make certain models more efficient by using ldpf directly
      if (dist == 6) { // lognormal
        for (i in 1:n) l[i] = lognormal_lpdf(time | eta[i], aux);
      } else if (dist == 8) { // Gamma
        vector[n] eeta = exp(-eta);
        for (i in 1:n) l[i] = gamma_lpdf(time | aux, eeta[i]);
      } else if (dist == 9) { // Gen Gamma
        for (i in 1:n) l[i] = gengamma_lpdf(time | eta[i], aux, aux2);
      } else {
        l = lS_a(dist, time, eta, aux, aux2) + lh_a(dist, time, eta, aux, aux2);
      }
    } else if (status == 2) { // Left censored
      l = log1m_exp(lS_a(dist, time, eta, aux, aux2));
    } else if (status == 3) { // Interval censored
      // l = log_diff_exp(lS_a(dist, start_time, eta, aux, aux2), lS_a(dist, time, eta, aux, aux2));
      l = log(exp(lS_a(dist, start_time, eta, aux, aux2)) - exp(lS_a(dist, time, eta, aux, aux2)));
    }

    // Left truncation
    if (delay_time > 0) {
      l -= lS_a(dist, delay_time, eta, aux, aux2);
    }

    return l;
  }
}
data {
#include /include/data_common.stan

  // Prior on aux parameter (usually shape)
  int<lower=0,upper=6> prior_aux_dist;
  real prior_aux_location;
  real<lower=0> prior_aux_scale;
  real<lower=0> prior_aux_df;

  // Prior on aux2 parameter for gengamma
  int<lower=0,upper=6> prior_aux2_dist;
  real prior_aux2_location;
  real<lower=0> prior_aux2_scale;
  real<lower=0> prior_aux2_df;

  // Select distribution
  int<lower=1> dist;
  // 1 = Exponential PH
  // 2 = Weibull PH
  // 3 = Gompertz PH
  // 4 = Exponential AFT
  // 5 = Weibull AFT
  // 6 = log Normal AFT
  // 7 = log logistic AFT
  // 8 = Gamma AFT
  // 9 = Generalised Gamma AFT

  // AgD arm-based
  // int<lower=0> narm_agd_arm; // Number of arm-based AgD arms
  int<lower=1> agd_arm_arm[ni_agd_arm]; // Arm indicator for AgD (arm-based) (i.e. picking element of which_RE)

  // Outcomes
  vector[ni_ipd] ipd_time;
  vector[ni_ipd] ipd_start_time;
  vector[ni_ipd] ipd_delay_time;
  int<lower=0, upper=3> ipd_status[ni_ipd];

  vector[ni_agd_arm] agd_arm_time;
  vector[ni_agd_arm] agd_arm_start_time;
  vector[ni_agd_arm] agd_arm_delay_time;
  int<lower=0, upper=3> agd_arm_status[ni_agd_arm];

  // Aux IDs for independent shape parameters
  int<lower=0, upper=1> aux_by; // Flag subgroup aux parameters within each arm (1 = yes)
  int<lower=1> aux_id[ni_ipd + ni_agd_arm * (aux_by ? nint_max : 1)];
}
transformed data {
  // Exponential model indicator, 0 = exponential
  int<lower=0, upper=1> nonexp = (dist == 1 || dist == 4) ? 0 : 1;
  // Generalised Gamma model indicator, 1 = gengamma
  int<lower=0, upper=1> gengamma = dist >= 9 ? 1 : 0;
  // Number of auxiliary parameters
  int n_aux = nonexp ? max(aux_id) : 0;
  // AgD arm indicator - shifted
  int agd_arm_arm2[ni_agd_arm];
  // Number of integration points

#include /include/transformed_data_common.stan

  for (i in 1:ni_agd_arm) agd_arm_arm2[i] = narm_ipd + agd_arm_arm[i];
}
parameters {
#include /include/parameters_common.stan

  // Auxiliary parameters for parametric model (typically shape)
  // Exponential model has shape = 1 so parameter is removed (zero dim)
  vector<lower=0>[n_aux] aux;
  // Second auxiliary parameter for generalised gamma
  vector<lower=0>[n_aux*gengamma] aux2;
}
transformed parameters {
  // Log likelihood contributions
  vector[ni_ipd] log_L_ipd;
  vector[ni_agd_arm] log_L_agd_arm;

#include /include/transformed_parameters_common.stan


  // Evaluate log likelihood
  log_L_ipd = loglik(dist,
                     ipd_time,
                     ipd_start_time,
                     ipd_delay_time,
                     ipd_status,
                     eta_ipd,
                     nonexp ? aux[aux_id[1:ni_ipd]] : rep_vector(0, ni_ipd),
                     gengamma ? aux2[aux_id[1:ni_ipd]] : rep_vector(0, ni_ipd));

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
        if (RE && which_RE[agd_arm_arm2[i]]) eta_agd_arm_ii += f_delta[which_RE[agd_arm_arm2[i]]];

        // Average likelihood over integration points
        if (aux_by) {
          log_L_ii = loglik(dist,
                            rep_vector(agd_arm_time[i], nint),
                            rep_vector(agd_arm_start_time[i], nint),
                            rep_vector(agd_arm_delay_time[i], nint),
                            rep_array(agd_arm_status[i], nint),
                            eta_agd_arm_ii,
                            nonexp ? aux[aux_id[(ni_ipd + 1 + (i-1)*nint_max):(ni_ipd + (i-1)*nint_max + nint)]] : rep_vector(0, nint),
                            gengamma ? aux2[aux_id[(ni_ipd + 1 + (i-1)*nint_max):(ni_ipd + (i-1)*nint_max + nint)]] : rep_vector(0, nint));
        } else {
          log_L_ii = loglik_a(dist,
                              agd_arm_time[i],
                              agd_arm_start_time[i],
                              agd_arm_delay_time[i],
                              agd_arm_status[i],
                              eta_agd_arm_ii,
                              nonexp ? aux[aux_id[ni_ipd + i]] : 0,
                              gengamma ? aux2[aux_id[ni_ipd + i]] : 0);
        }

        log_L_agd_arm[i] = log_sum_exp(log_L_ii) - log(nint);
      }
    } else { // -- If no integration --
      vector[ni_agd_arm] eta_agd_arm = eta_agd_arm_noRE;

      // Random effects
      if (RE) for (i in 1:ni_agd_arm) if (which_RE[agd_arm_arm2[i]]) eta_agd_arm[i] += f_delta[which_RE[agd_arm_arm2[i]]];

      log_L_agd_arm = loglik(dist,
                             agd_arm_time,
                             agd_arm_start_time,
                             agd_arm_delay_time,
                             agd_arm_status,
                             eta_agd_arm,
                             nonexp ? aux[aux_id[(ni_ipd + 1):(ni_ipd + ni_agd_arm)]] : rep_vector(0, ni_agd_arm),
                             gengamma ? aux2[aux_id[(ni_ipd + 1):(ni_ipd + ni_agd_arm)]] : rep_vector(0, ni_agd_arm));
    }
  }
}
model {
#include /include/model_common.stan

  // -- Prior on auxiliary parameters --
  if (nonexp) prior_select_lp(aux, prior_aux_dist, prior_aux_location, prior_aux_scale, prior_aux_df);
  if (gengamma) prior_select_lp(aux2, prior_aux2_dist, prior_aux2_location, prior_aux2_scale, prior_aux2_df);

  // -- IPD likelihood --
  target += log_L_ipd;

  // -- AgD likelihood (arm-based) --
  target += log_L_agd_arm;

}
generated quantities {
  // Transform intercepts back to scales
  // vector[ns_ipd + ns_agd_arm] scale = exp(mu);

  vector[(dist != 1 && dist != 4 && dist != 6 && dist != 9) ? n_aux : 0] shape;
  vector[dist == 6 ? n_aux : 0] sdlog;  // lognormal sdlog
  vector[dist == 9 ? n_aux : 0] sigma; // gengamma sigma
  vector[dist == 9 ? n_aux : 0] k;  // gengamma k

#include /include/generated_quantities_common.stan

  // Log likelihood
  if (ni_ipd) log_lik[1:ni_ipd] = log_L_ipd;
  if (ni_agd_arm) log_lik[(ni_ipd + 1):(ni_ipd + ni_agd_arm)] = log_L_agd_arm;

  // Residual deviance
  // NOTE: This is actually just the deviance, which is fine for relative fit model comparison
  if (ni_ipd + ni_agd_arm) resdev[1:(ni_ipd + ni_agd_arm)] = -2 * log_lik[1:(ni_ipd + ni_agd_arm)];

  // Fitted values not implemented

  // Rename parameters
  if (dist != 1 && dist != 4 && dist != 6 && dist != 9) shape = aux;
  if (dist == 6) sdlog = aux;
  if (dist == 9) {
    sigma = aux;
    k = aux2;
  }

}
