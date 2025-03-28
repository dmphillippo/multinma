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
      out = -pow(y, aux) .* exp(eta);
    } else if (dist == 3) { // Gompertz
      out = -exp(eta)./aux .* expm1(aux .* y);
    } else if (dist == 4) { // Exponential AFT
      out = -y .* exp(-eta);
    } else if (dist == 5) { // Weibull AFT
      out = -pow(y, aux) .* exp(-aux .* eta);
    } else if (dist == 6) { // log Normal
      // aux is sdlog
      // out = log1m(Phi((log(y) - eta) ./ aux));
      for (i in 1:n) out[i] = lognormal_lccdf(y[i] | eta[i], aux[i]);
    } else if (dist == 7) { // log logistic
      out = -log1p(pow(y ./ exp(eta), aux));
    } else if (dist == 8) { // Gamma
      vector[n] eeta = exp(-eta);
      for (i in 1:n) out[i] = gamma_lccdf(y[i] | aux[i], eeta[i]);
    } else if (dist == 9) { // Generalised Gamma
      vector[n] Q = inv(sqrt(aux2));
      vector[n] w = exp(Q .* (log(y) - eta) ./ aux) .* aux2;
      out = log1m(gamma_p(aux2, w));
    }

    return out;
  }

  // Version with single aux parameter for all individuals
  vector lS2(int dist, vector y, vector eta, real aux, real aux2) {
    // aux is shape, except for lognormal sdlog, gengamma sigma
    // aux2 is only for gengamma k

    int n = num_elements(y);
    vector[n] out;

    if (dist == 1) { // Exponential
      out = -y .* exp(eta);
    } else if (dist == 2) { // Weibull
      out = -pow(y, aux) .* exp(eta);
    } else if (dist == 3) { // Gompertz
      out = -exp(eta)/aux .* expm1(aux * y);
    } else if (dist == 4) { // Exponential AFT
      out = -y .* exp(-eta);
    } else if (dist == 5) { // Weibull AFT
      out = -pow(y, aux) .* exp(-aux * eta);
    } else if (dist == 6) { // log Normal
      // aux is sdlog
      // out = log1m(Phi((log(y) - eta) ./ aux));
      for (i in 1:n) out[i] = lognormal_lccdf(y[i] | eta[i], aux);
    } else if (dist == 7) { // log logistic
      out = -log1p(pow(y ./ exp(eta), aux));
    } else if (dist == 8) { // Gamma
      vector[n] eeta = exp(-eta);
      for (i in 1:n) out[i] = gamma_lccdf(y[i] | aux, eeta[i]);
    } else if (dist == 9) { // Generalised Gamma
      real Q = inv(sqrt(aux2));
      vector[n] w = exp(Q * (log(y) - eta) / aux) * aux2;
      out = log1m(gamma_p(aux2, w));
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
      out = -log1p(pow(y ./ exp(eta), aux));
    } else if (dist == 8) { // Gamma
      vector[n] eeta = exp(-eta);
      for (i in 1:n) out[i] = gamma_lccdf(y | aux, eeta[i]);
    } else if (dist == 9) { // Generalised Gamma
      real Q = pow(aux2, -0.5);
      vector[n] w = exp(Q * (log(y) - eta) / aux) * aux2;
      out = log1m(gamma_p(aux2, w));
    }

    return out;
  }

  // AgD version with integration points over eta and aux pars
  vector lS_a2(int dist, real y, vector eta, vector aux, vector aux2) {
    // aux is shape, except for lognormal sdlog, gengamma sigma
    // aux2 is only for gengamma k

    int n = num_elements(eta);
    vector[n] out;

    if (dist == 1) { // Exponential
      out = -y * exp(eta);
    } else if (dist == 2) { // Weibull
      out = -pow(y, aux) .* exp(eta);
    } else if (dist == 3) { // Gompertz
      out = -exp(eta)./aux .* expm1(aux * y);
    } else if (dist == 4) { // Exponential AFT
      out = -y * exp(-eta);
    } else if (dist == 5) { // Weibull AFT
      out = -pow(y, aux) .* exp(-aux .* eta);
    } else if (dist == 6) { // log Normal
      // aux is sdlog
      // lS = log1m(Phi((log(y) - eta) / aux));
      for (i in 1:n) out[i] = lognormal_lccdf(y | eta[i], aux[i]);
    } else if (dist == 7) { // log logistic
      out = -log1p(pow(y ./ exp(eta), aux));
    } else if (dist == 8) { // Gamma
      vector[n] eeta = exp(-eta);
      for (i in 1:n) out[i] = gamma_lccdf(y | aux[i], eeta[i]);
    } else if (dist == 9) { // Generalised Gamma
      vector[n] Q = pow(aux2, -0.5);
      vector[n] w = exp(Q .* (log(y) - eta) ./ aux) .* aux2;
      out = log1m(gamma_p(aux2, w));
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
      out = log(aux) + eta + lmultiply(aux - 1, y);
    } else if (dist == 3) { // Gompertz
      out = eta + (aux .* y);
    } else if (dist == 4) { // Exponential AFT
      out = -eta;
    } else if (dist == 5) { // Weibull AFT
      out = log(aux) - (aux .* eta) + lmultiply(aux - 1, y);
    } else if (dist == 6) { // log Normal
      // aux is sdlog
      // out = lognormal_lpdf(y | eta, aux) - log1m(Phi((log(y) - eta) ./ aux));
      for (i in 1:n) out[i] = lognormal_lpdf(y[i] | eta[i], aux[i]) - lognormal_lccdf(y[i] | eta[i], aux[i]);
    } else if (dist == 7) { // log logistic
      out = log(aux) - eta + (aux - 1).*(log(y) - eta) - log1p(pow(y ./ exp(eta), aux));
    } else if (dist == 8) { // Gamma
      vector[n] eeta = exp(-eta);
      for (i in 1:n) out[i] = gamma_lpdf(y[i] | aux[i], eeta[i]) - gamma_lccdf(y[i] | aux[i], eeta[i]);
    } else if (dist == 9) { // Generalised Gamma
      // Not used, lpdf used directly
    }

    return out;
  }

  // Version with single aux parameter for all individuals
  vector lh2(int dist, vector y, vector eta, real aux, real aux2) {
    // aux is shape, except for lognormal sdlog, gengamma scale
    // aux2 is only for gengamma shape

    int n = num_elements(y);
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
      for (i in 1:n) out[i] = lognormal_lpdf(y[i] | eta[i], aux) - lognormal_lccdf(y[i] | eta[i], aux);
    } else if (dist == 7) { // log logistic
      out = log(aux) - eta + (aux - 1)*(log(y) - eta) - log1p(pow(y ./ exp(eta), aux));
    } else if (dist == 8) { // Gamma
      vector[n] eeta = exp(-eta);
      for (i in 1:n) out[i] = gamma_lpdf(y[i] | aux, eeta[i]) - gamma_lccdf(y[i] | aux, eeta[i]);
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
      out = log(aux) - eta + (aux - 1)*(log(y) - eta) - log1p(pow(y ./ exp(eta), aux));
    } else if (dist == 8) { // Gamma
      vector[n] eeta = exp(-eta);
      for (i in 1:n) out[i] = gamma_lpdf(y | aux, eeta[i]) - gamma_lccdf(y | aux, eeta[i]);
    } else if (dist == 9) { // Generalised Gamma
      // Not used, lpdf used directly
    }

    return out;
  }

  // AgD version with integration points over eta and aux pars
  vector lh_a2(int dist, real y, vector eta, vector aux, vector aux2) {
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
      out = log(aux) - (aux .* eta) + lmultiply(aux - 1, y);
    } else if (dist == 6) { // log Normal
      // aux is sdlog
      // out = lognormal_lpdf(y | eta, aux) - log1m(Phi((log(y) - eta) ./ aux));
      for (i in 1:n) out[i] = lognormal_lpdf(y| eta[i], aux[i]) - lognormal_lccdf(y | eta[i], aux[i]);
    } else if (dist == 7) { // log logistic
      out = log(aux) - eta + (aux - 1).*(log(y) - eta) - log1p(pow(y ./ exp(eta), aux));
    } else if (dist == 8) { // Gamma
      vector[n] eeta = exp(-eta);
      for (i in 1:n) out[i] = gamma_lpdf(y | aux[i], eeta[i]) - gamma_lccdf(y | aux[i], eeta[i]);
    } else if (dist == 9) { // Generalised Gamma
      // Not used, lpdf used directly
    }

    return out;
  }

  // -- Log likelihood with censoring and truncation --
  vector loglik(int dist, vector time, vector start_time, vector delay_time, array[] int status, vector eta, vector aux, vector aux2) {
    int n = num_elements(eta);
    vector[n] l;

    // Status indicators
    array[3] int nw = nwhich_all(status, 3);
    int nw0 = (dist == 6 || dist == 8 || dist == 9) ? n - sum(nw) : 0;
    array[nw0] int w0;
    array[nw[1]] int w1;
    array[nw[2]] int w2;
    array[nw[3]] int w3;
    int nd = num_elements(which_gt0(delay_time));
    array[nd] int wd;
    if (nw0) w0 = which(status, 0);
    if (nw[1]) w1 = which(status, 1);
    if (nw[2]) w2 = which(status, 2);
    if (nw[3]) w3 = which(status, 3);
    if (nd) wd = which_gt0(delay_time);

    // Make certain models more efficient by using ldpf directly
    if (dist == 6 || dist == 8 || dist == 9) {
      // Right censored
      l[w0] = lS(dist, time[w0], eta[w0], aux[w0], aux2[w0]);

      // Observed
      if (dist == 6) for (i in 1:n) if (status[i] == 1) l[i] = lognormal_lpdf(time[i] | eta[i], aux[i]);
      if (dist == 8) for (i in 1:n) if (status[i] == 1) l[i] = gamma_lpdf(time[i] | aux[i], exp(-eta[i]));
      if (dist == 9) for (i in 1:n) if (status[i] == 1) l[i] = gengamma_lpdf(time[i] | eta[i], aux[i], aux2[i]);

    } else {

      // Right censored
      l = lS(dist, time, eta, aux, aux2);

      // Observed
      l[w1] += lh(dist, time[w1], eta[w1], aux[w1], aux2[w1]);
    }

    // Left censored
    l[w2] = log1m_exp(l[w2]);

    // Interval censored
    // l = log_diff_exp(lS(dist, start_time, eta, aux, aux2), lS(dist, time, eta, aux, aux2));
    // l[w3] = log(exp(lS(dist, start_time[w3], eta[w3], aux[w3], aux2[w3])) - exp(l[w3]));
    l[w3] = log_diff_exp(lS(dist, start_time[w3], eta[w3], aux[w3], aux2[w3]), l[w3]);

    // Left truncation
    l[wd] -= lS(dist, delay_time[wd], eta[wd], aux[wd], aux2[wd]);

    return l;
  }

  // Version with single aux parameter for all individuals
  vector loglik2(int dist, vector time, vector start_time, vector delay_time, array[] int status, vector eta, real aux, real aux2) {
    int n = num_elements(eta);
    vector[n] l;

    // Status indicators
    array[3] int nw = nwhich_all(status, 3);
    int nw0 = (dist == 6 || dist == 8 || dist == 9) ? n - sum(nw) : 0;
    array[nw0] int w0;
    array[nw[1]] int w1;
    array[nw[2]] int w2;
    array[nw[3]] int w3;
    int nd = num_elements(which_gt0(delay_time));
    array[nd] int wd;
    if (nw0) w0 = which(status, 0);
    if (nw[1]) w1 = which(status, 1);
    if (nw[2]) w2 = which(status, 2);
    if (nw[3]) w3 = which(status, 3);
    if (nd) wd = which_gt0(delay_time);

    // Make certain models more efficient by using ldpf directly
    if (dist == 6 || dist == 8 || dist == 9) {
      // Right censored
      l[w0] = lS2(dist, time[w0], eta[w0], aux, aux2);

      // Observed
      if (dist == 6) for (i in 1:n) if (status[i] == 1) l[i] = lognormal_lpdf(time[i] | eta[i], aux);
      if (dist == 8) for (i in 1:n) if (status[i] == 1) l[i] = gamma_lpdf(time[i] | aux, exp(-eta[i]));
      if (dist == 9) for (i in 1:n) if (status[i] == 1) l[i] = gengamma_lpdf(time[i] | eta[i], aux, aux2);

    } else {

      // Right censored
      l = lS2(dist, time, eta, aux, aux2);

      // Observed
      l[w1] += lh2(dist, time[w1], eta[w1], aux, aux2);
    }

    // Left censored
    l[w2] = log1m_exp(l[w2]);

    // Interval censored
    // l = log_diff_exp(lS(dist, start_time, eta, aux, aux2), lS(dist, time, eta, aux, aux2));
    // l[w3] = log(exp(lS(dist, start_time[w3], eta[w3], aux[w3], aux2[w3])) - exp(l[w3]));
    l[w3] = log_diff_exp(lS2(dist, start_time[w3], eta[w3], aux, aux2), l[w3]);

    // Left truncation
    l[wd] -= lS2(dist, delay_time[wd], eta[wd], aux, aux2);

    return l;
  }

  // AgD version with integration points over eta
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
      l = log_diff_exp(lS_a(dist, start_time, eta, aux, aux2), lS_a(dist, time, eta, aux, aux2));
      // l = log(exp(lS_a(dist, start_time, eta, aux, aux2)) - exp(lS_a(dist, time, eta, aux, aux2)));
    }

    // Left truncation
    if (delay_time > 0) {
      l -= lS_a(dist, delay_time, eta, aux, aux2);
    }

    return l;
  }

  // AgD version with integration points over eta and aux pars
  vector loglik_a2(int dist, real time, real start_time, real delay_time, int status, vector eta, vector aux, vector aux2) {
    int n = num_elements(eta);
    vector[n] l;

    if (status == 0) { // Right censored
      l = lS_a2(dist, time, eta, aux, aux2);
    } else if (status == 1) { // Observed
      // Make certain models more efficient by using ldpf directly
      if (dist == 6) { // lognormal
        for (i in 1:n) l[i] = lognormal_lpdf(time | eta[i], aux[i]);
      } else if (dist == 8) { // Gamma
        vector[n] eeta = exp(-eta);
        for (i in 1:n) l[i] = gamma_lpdf(time | aux[i], eeta[i]);
      } else if (dist == 9) { // Gen Gamma
        for (i in 1:n) l[i] = gengamma_lpdf(time | eta[i], aux[i], aux2[i]);
      } else {
        l = lS_a2(dist, time, eta, aux, aux2) + lh_a2(dist, time, eta, aux, aux2);
      }
    } else if (status == 2) { // Left censored
      l = log1m_exp(lS_a2(dist, time, eta, aux, aux2));
    } else if (status == 3) { // Interval censored
      l = log_diff_exp(lS_a2(dist, start_time, eta, aux, aux2), lS_a2(dist, time, eta, aux, aux2));
      // l = log(exp(lS_a2(dist, start_time, eta, aux, aux2)) - exp(lS_a2(dist, time, eta, aux, aux2)));
    }

    // Left truncation
    if (delay_time > 0) {
      l -= lS_a2(dist, delay_time, eta, aux, aux2);
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
  array[ni_agd_arm] int<lower=1> agd_arm_arm; // Arm indicator for AgD (arm-based) (i.e. picking element of which_RE)

  // Outcomes
  vector[ni_ipd] ipd_time;
  vector[ni_ipd] ipd_start_time;
  vector[ni_ipd] ipd_delay_time;
  array[ni_ipd] int<lower=0, upper=3> ipd_status;

  vector[ni_agd_arm] agd_arm_time;
  vector[ni_agd_arm] agd_arm_start_time;
  vector[ni_agd_arm] agd_arm_delay_time;
  array[ni_agd_arm] int<lower=0, upper=3> agd_arm_status;

  // Aux IDs for independent shape parameters
  int<lower=0, upper=1> aux_int; // Flag aux regression needs integration (1 = yes)
  //int<lower=0, upper=1> aux_by; // Flag subgroup aux parameters within each arm (1 = yes)
  array[ni_ipd + ni_agd_arm * (aux_int ? nint_max : 1)] int<lower=1> aux_id;
  array[ni_ipd + ni_agd_arm * (aux_int ? nint_max : 1)] int<lower=1> aux_group;

  // auxiliary design matrix
  int<lower=0> nX_aux;
  matrix[ni_ipd + (aux_int ? nint_max : 1) * ni_agd_arm, nX_aux] X_aux;
  int<lower=0, upper=1> aux_reg_trt; // Flag aux regression includes main treatment effects (1 = yes)

  // Prior on aux regression coefficients
  int<lower=0,upper=6> prior_aux_reg_dist;
  real<lower=0> prior_aux_reg_location;
  real<lower=0> prior_aux_reg_scale;
  real<lower=0> prior_aux_reg_df;
}
transformed data {
  // Exponential model indicator, 0 = exponential
  int<lower=0, upper=1> nonexp = (dist == 1 || dist == 4) ? 0 : 1;
  // Generalised Gamma model indicator, 1 = gengamma
  int<lower=0, upper=1> gengamma = dist >= 9 ? 1 : 0;
  // Number of auxiliary parameters
  int n_aux = nonexp ? max(aux_id) : 0;
  int n_aux_group = nonexp ? max(aux_group) : 0;
  // AgD arm indicator - shifted
  array[ni_agd_arm] int agd_arm_arm2;

  // Aux IDs
  array[ni_ipd] int<lower=1> aux_id_ipd = aux_id[1:ni_ipd];
  array[ni_agd_arm*(aux_int ? nint_max : 1)] int<lower=1> aux_id_agd_arm = aux_id[(ni_ipd + 1):(ni_ipd + ni_agd_arm*(aux_int ? nint_max : 1))];

  array[ni_ipd] int<lower=1> aux_group_ipd = aux_group[1:ni_ipd];
  array[ni_agd_arm*(aux_int ? nint_max : 1)] int<lower=1> aux_group_agd_arm = aux_group[(ni_ipd + 1):(ni_ipd + ni_agd_arm*(aux_int ? nint_max : 1))];

  // Aux ID indexing arrays
  array[(nonexp && aux_int == 0) ? n_aux_group : 1] int ni_aux_group_ipd = (nonexp && aux_int == 0) ? nwhich_all(aux_group_ipd, n_aux_group) : {0};
  array[(nonexp && aux_int == 0) ? n_aux_group : 0, (nonexp && aux_int == 0) ? max(ni_aux_group_ipd) : 0] int wi_aux_group_ipd;
  array[(nonexp && aux_int == 0) ? n_aux_group : 1] int ni_aux_group_agd_arm = (nonexp && aux_int == 0) ? nwhich_all(aux_group_agd_arm, n_aux_group) : {0};
  array[(nonexp && aux_int == 0) ? n_aux_group : 0, (nonexp && aux_int == 0) ? max(ni_aux_group_agd_arm) : 0] int wi_aux_group_agd_arm;

  // Split auxiliary design matrix into IPD and AgD rows
  matrix[0, nX_aux] Xauxdummy;
  matrix[ni_ipd, nX_aux] X_aux_ipd = ni_ipd ? X_aux[1:ni_ipd] : Xauxdummy;
  matrix[(aux_int ? nint_max : 1) * ni_agd_arm, nX_aux] X_aux_agd_arm = ni_agd_arm ? X_aux[(ni_ipd + 1):(ni_ipd + (aux_int ? nint_max : 1) * ni_agd_arm)] : Xauxdummy;

#include /include/transformed_data_common.stan

  for (i in 1:ni_agd_arm) agd_arm_arm2[i] = narm_ipd + agd_arm_arm[i];

  if (nonexp && aux_int == 0) for (i in 1:n_aux_group) {
    if (ni_aux_group_ipd[i]) wi_aux_group_ipd[i, 1:ni_aux_group_ipd[i]] = which(aux_group_ipd, i);
    if (ni_aux_group_agd_arm[i]) wi_aux_group_agd_arm[i, 1:ni_aux_group_agd_arm[i]] = which(aux_group_agd_arm, i);
  }
}
parameters {
#include /include/parameters_common.stan

  // Auxiliary parameters for parametric model (typically shape)
  // Exponential model has shape = 1 so parameter is removed (zero dim)
  vector<lower=0>[n_aux] aux;
  // Second auxiliary parameter for generalised gamma
  vector<lower=0>[n_aux*gengamma] aux2;

  // Auxiliary regression
  matrix[nonexp ? nX_aux : 0, nonexp + gengamma] beta_aux;
}
transformed parameters {
  // Log likelihood contributions
  vector[ni_ipd] log_L_ipd;
  vector[ni_agd_arm] log_L_agd_arm;

#include /include/transformed_parameters_common.stan


  // Evaluate log likelihood
  if (ni_ipd) {
    if (nonexp == 0 || aux_int) {
      vector[ni_ipd] auxi;
      vector[ni_ipd] aux2i;
      matrix[nonexp ? ni_ipd : 0, nonexp + gengamma] eXbeta;
      if (nonexp) eXbeta = exp(X_aux_ipd * beta_aux);
      if (nonexp) auxi = aux[aux_id_ipd] .* eXbeta[,1];
      if (gengamma) aux2i = aux2[aux_id_ipd] .* eXbeta[,2];

      log_L_ipd = loglik(dist,
                         ipd_time,
                         ipd_start_time,
                         ipd_delay_time,
                         ipd_status,
                         eta_ipd,
                         auxi,
                         aux2i);
    } else {
      for (i in 1:n_aux_group) {
        int ni = ni_aux_group_ipd[i];

        if (ni) {
          array[ni] int wi = wi_aux_group_ipd[i, 1:ni];
          real auxi;
          real aux2i;
          row_vector[nonexp + gengamma] eXbeta;
          if (nonexp) eXbeta = exp(X_aux_ipd[wi[1],] * beta_aux);

          if (nX_aux) {
            if (nonexp) auxi = aux[aux_id_ipd[wi[1]]] * eXbeta[1];
            if (gengamma) aux2i = aux2[aux_id_ipd[wi[1]]] * eXbeta[2];
          } else {
            // If no X_aux then aux_id = aux_group
            if (nonexp) auxi = aux[i];
            if (gengamma) aux2i = aux2[i];
          }

          log_L_ipd[wi] = loglik2(dist,
                                  ipd_time[wi],
                                  ipd_start_time[wi],
                                  ipd_delay_time[wi],
                                  ipd_status[wi],
                                  eta_ipd[wi],
                                  auxi,
                                  aux2i);
        }
      }
    }
  }

  // -- AgD model (arm-based) --
  if (ni_agd_arm) {
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

    if (nint_max > 1) { // -- If integration points are used --

      if (nonexp == 0 || aux_int) {
        for (i in 1:ni_agd_arm) {
          vector[nint] eta_agd_arm_ii;
          vector[nint] log_L_ii;
          vector[nint] auxi;
          vector[nint] aux2i;
          matrix[nonexp ? nint : 0, nonexp + gengamma] eXbeta;

          if (nonexp) eXbeta = exp(X_aux_agd_arm[(1 + (i-1)*nint_max):((i-1)*nint_max + nint), ] * beta_aux);
          if (nonexp) auxi = aux[aux_id_agd_arm[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)]] .* eXbeta[,1];
          if (gengamma) aux2i = aux2[aux_id_agd_arm[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)]] .* eXbeta[,2];

          eta_agd_arm_ii = eta_agd_arm_noRE[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)];

          // Random effects
          if (RE && which_RE[agd_arm_arm2[i]]) eta_agd_arm_ii += f_delta[which_RE[agd_arm_arm2[i]]];

          // Average likelihood over integration points
          log_L_ii = loglik_a2(dist,
                               agd_arm_time[i],
                               agd_arm_start_time[i],
                               agd_arm_delay_time[i],
                               agd_arm_status[i],
                               eta_agd_arm_ii,
                               auxi,
                               aux2i);


          log_L_agd_arm[i] = log_sum_exp(log_L_ii) - log(nint);
        }
      } else {
        for (i in 1:n_aux_group) {
          int ni = ni_aux_group_agd_arm[i];

          if (ni) {
            array[ni] int wi = wi_aux_group_agd_arm[i, 1:ni];
            real auxi;
            real aux2i;
            row_vector[nonexp + gengamma] eXbeta;

            if (nonexp) eXbeta = exp(X_aux_agd_arm[wi[1],] * beta_aux);

            if (nX_aux) {
              auxi = nonexp ? aux[aux_id_agd_arm[wi[1]]] * eXbeta[1] : 0;
              aux2i = gengamma ? aux2[aux_id_agd_arm[wi[1]]] * eXbeta[2] : 0;
            } else {
              // If no X_aux then aux_id = aux_group
              auxi = nonexp ? aux[i] : 0;
              aux2i = gengamma ? aux2[i] : 0;
            }

            for (j in 1:ni) {
              vector[nint] eta_agd_arm_ii;
              vector[nint] log_L_ii;

              eta_agd_arm_ii = eta_agd_arm_noRE[(1 + (wi[j]-1)*nint_max):((wi[j]-1)*nint_max + nint)];
              if (RE && which_RE[agd_arm_arm2[wi[j]]]) eta_agd_arm_ii += f_delta[which_RE[agd_arm_arm2[wi[j]]]];

              log_L_ii = loglik_a(dist,
                                  agd_arm_time[wi[j]],
                                  agd_arm_start_time[wi[j]],
                                  agd_arm_delay_time[wi[j]],
                                  agd_arm_status[wi[j]],
                                  eta_agd_arm_ii,
                                  auxi,
                                  aux2i);

              log_L_agd_arm[wi[j]] = log_sum_exp(log_L_ii) - log(nint);
            }
          }
        }
      }
    } else { // -- If no integration --
      vector[ni_agd_arm] eta_agd_arm = eta_agd_arm_noRE;

      // Random effects
      if (RE) for (i in 1:ni_agd_arm) if (which_RE[agd_arm_arm2[i]]) eta_agd_arm[i] += f_delta[which_RE[agd_arm_arm2[i]]];

      if (nonexp == 0 || aux_int) {
        vector[ni_agd_arm] auxi;
        vector[ni_agd_arm] aux2i;
        matrix[ni_agd_arm, nonexp + gengamma] eXbeta;
        if (nonexp) eXbeta = exp(X_aux_agd_arm * beta_aux);
        if (nonexp) auxi = aux[aux_id_agd_arm] .* eXbeta[,1];
        if (gengamma) aux2i = aux2[aux_id_agd_arm] .* eXbeta[,2];

        log_L_agd_arm = loglik(dist,
                               agd_arm_time,
                               agd_arm_start_time,
                               agd_arm_delay_time,
                               agd_arm_status,
                               eta_agd_arm,
                               auxi,
                               aux2i);
      } else {
        for (i in 1:n_aux_group) {
          int ni = ni_aux_group_agd_arm[i];

          if (ni) {
            array[ni] int wi = wi_aux_group_agd_arm[i, 1:ni];
            real auxi;
            real aux2i;
            row_vector[nonexp + gengamma] eXbeta;

            if (nonexp) eXbeta = exp(X_aux_agd_arm[wi[1],] * beta_aux);

            if (nX_aux) {
              if (nonexp) auxi = aux[aux_id_agd_arm[wi[1]]] * eXbeta[1];
              if (gengamma) aux2i = aux2[aux_id_agd_arm[wi[1]]] * eXbeta[2];
            } else {
              // If no X_aux then aux_id = aux_group
              if (nonexp) auxi = aux[i];
              if (gengamma) aux2i = aux2[i];
            }

            log_L_agd_arm[wi] = loglik2(dist,
                                        agd_arm_time[wi],
                                        agd_arm_start_time[wi],
                                        agd_arm_delay_time[wi],
                                        agd_arm_status[wi],
                                        eta_agd_arm[wi],
                                        auxi,
                                        aux2i);
          }
        }
      }
    }
  }
}
model {
#include /include/model_common.stan

  // -- Prior on auxiliary parameters --
  if (nonexp) prior_select_lp(aux, prior_aux_dist, prior_aux_location, prior_aux_scale, prior_aux_df);
  if (gengamma) prior_select_lp(aux2, prior_aux2_dist, prior_aux2_location, prior_aux2_scale, prior_aux2_df);

  // -- Prior on aux regression --
  if (nonexp && nX_aux) prior_select_lp(beta_aux[,1], prior_aux_reg_dist, prior_aux_reg_location, prior_aux_reg_scale, prior_aux_reg_df);
  if (gengamma && nX_aux) prior_select_lp(beta_aux[,2], prior_aux_reg_dist, prior_aux_reg_location, prior_aux_reg_scale, prior_aux_reg_df);

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

  // aux regression treatment effects
  matrix[aux_reg_trt ? nt-1 : 0,  nonexp + gengamma] d_aux;

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

  // aux regression treatment effects
  if (aux_reg_trt) d_aux = beta_aux[1:(nt-1), ];

}
