functions {
#include /include/prior_select.stan
#include /include/count_nonzero.stan
#include /include/which.stan

  // -- Vector functions --
  vector pow_vec (vector x, vector y) {
    // (Not needed past Stan 2.26)
    int n = num_elements(x);
    vector[n] out;
    for (i in 1:n) out[i] = pow(x[i], y[i]);
    return out;
  }

  vector lmultiply_vec (vector x, vector y) {
    int n = num_elements(x);
    vector[n] out;
    for (i in 1:n) out[i] = lmultiply(x[i], y[i]);
    return out;
  }

  // Which entries are greater than 0
  int[] which_gt0(vector x) {
    int n = num_elements(x);
    int w[n]; // Over-allocate w and then truncate later
    int c = 1;
    for (i in 1:n) {
      if (x[i] > 0) {
        w[c] = i;
        c += 1;
      }
    }
    return w[1:(c-1)];
  }

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
    vector[n] lS;

    if (dist == 1) { // Exponential
      lS = -y .* exp(eta);
    } else if (dist == 2) { // Weibull
      lS = -pow_vec(y, aux) .* exp(eta);
    } else if (dist == 3) { // Gompertz
      lS = -exp(eta)./aux .* expm1(aux .* y);
    } else if (dist == 4) { // Exponential AFT
      lS = -y .* exp(-eta);
    } else if (dist == 5) { // Weibull AFT
      lS = -pow_vec(y, aux) .* exp(-aux .* eta);
    } else if (dist == 6) { // log Normal
      // aux is sdlog
      // lS = log1m(Phi((log(y) - eta) ./ aux));
      for (i in 1:n) lS[i] = lognormal_lccdf(y[i] | eta[i], aux[i]);
    } else if (dist == 7) { // log logistic
      lS = -log1p(pow_vec(y ./ exp(eta), aux));
    } else if (dist == 8) { // Gamma
      for (i in 1:n) lS[i] = gamma_lccdf(y[i] | aux[i], exp(-eta[i]));
    } else if (dist == 9) { // Generalised Gamma
      vector[n] Q = inv(sqrt(aux2));
      vector[n] w = exp(Q .* (log(y) - eta) ./ aux) .* aux2;
      for (i in 1:n) lS[i] = log1m(gamma_p(aux2[i], w[i]));
    }

    return lS;
  }

  // -- log Hazard functions --
  vector lh(int dist, vector y, vector eta, vector aux, vector aux2) {
    // aux is shape, except for lognormal sdlog, gengamma scale
    // aux2 is only for gengamma shape

    int n = num_elements(y);
    vector[n] lh;

    if (dist == 1) { // Exponential
      lh = eta;
    } else if (dist == 2) { // Weibull
      lh = log(aux) + eta + lmultiply_vec(aux - 1, y);
    } else if (dist == 3) { // Gompertz
      lh = eta + (aux .* y);
    } else if (dist == 4) { // Exponential AFT
      lh = -eta;
    } else if (dist == 5) { // Weibull AFT
      lh = log(aux) - (aux .* eta) + lmultiply_vec(aux - 1, y);
    } else if (dist == 6) { // log Normal
      // aux is sdlog
      // lh = lognormal_lpdf(y | eta, aux) - log1m(Phi((log(y) - eta) ./ aux));
      for (i in 1:n) lh[i] = lognormal_lpdf(y[i] | eta[i], aux[i]) - lognormal_lccdf(y[i] | eta[i], aux[i]);
    } else if (dist == 7) { // log logistic
      lh = log(aux) - eta + (aux - 1).*(log(y) - eta) - log1p(pow_vec(y ./ exp(eta), aux));
    } else if (dist == 8) { // Gamma
      for (i in 1:n) lh[i] = gamma_lpdf(y[i] | aux[i], exp(-eta[i])) - gamma_lccdf(y[i] | aux[i], exp(-eta[i]));
    } else if (dist == 9) { // Generalised Gamma
      // Not used, lpdf used directly
    }

    return lh;
  }

  // -- Log likelihood with censoring and truncation --
  vector loglik(int dist, vector time, vector start_time, vector delay_time, int[] status, vector eta, vector aux, vector aux2) {
    int n = num_elements(eta);
    vector[n] l;

//     if (status == 0) { // Right censored
//       l = lS(dist, time, eta, aux, aux2);
//     } else if (status == 1) { // Observed
//       // Make (generalised) Gamma models more efficient by using ldpf directly
//       if (dist == 8) { // Gamma
//         l = gamma_lpdf(time | aux, exp(-eta));
//       } if (dist == 9) { // Gen Gamma
//         l = gengamma_lpdf(time | eta, aux, aux2);
//       } else {
//         l = lS(dist, time, eta, aux, aux2) + lh(dist, time, eta, aux, aux2);
//       }
//     } else if (status == 2) { // Left censored
//       l = log1m_exp(lS(dist, time, eta, aux, aux2));
//     } else if (status == 3) { // Interval censored
//       l = log_diff_exp(lS(dist, start_time, eta, aux, aux2), lS(dist, time, eta, aux, aux2));
//     }
//
//     // Left truncation
//     if (delay_time > 0) {
//       l -= lS(dist, delay_time, eta, aux, aux2);
//     }

    // Make certain models more efficient by using ldpf directly
    if (dist == 6 || dist == 8 || dist == 9) {
      // Right censored
      {
        int nwhich0 = num_elements(which(status, 0));
        int which0[nwhich0] = which(status, 0);
        l[which0] = lS(dist, time[which0], eta[which0], aux[which0], aux2[which0]);
      }

      // Observed
      if (dist == 6) for (i in 1:n) if (status[i] == 1) l[i] = lognormal_lpdf(time[i] | eta[i], aux[i]);
      if (dist == 8) for (i in 1:n) if (status[i] == 1) l[i] = gamma_lpdf(time[i] | aux[i], exp(-eta[i]));
      if (dist == 9) for (i in 1:n) if (status[i] == 1) l[i] = gengamma_lpdf(time[i] | eta[i], aux[i], aux2[i]);

    } else {

      // Right censored
      l = lS(dist, time, eta, aux, aux2);

      // Observed
      {
        int nwhich1 = num_elements(which(status, 1));
        int which1[nwhich1] = which(status, 1);
        l[which1] += lh(dist, time[which1], eta[which1], aux[which1], aux2[which1]);
      }
    }

    // Left censored
    l[which(status, 2)] = log1m_exp(l[which(status, 2)]);

    // Interval censored
    // l = log_diff_exp(lS(dist, start_time, eta, aux, aux2), lS(dist, time, eta, aux, aux2));
    {
      int nwhich3 = num_elements(which(status, 3));
      int which3[nwhich3] = which(status, 3);
      l[which3] = log(exp(l[which3]) - exp(lS(dist, time[which3], eta[which3], aux[which3], aux2[which3])));
    }

    // Left truncation
    {
      int nwhichd = num_elements(which_gt0(delay_time));
      int whichd[nwhichd] = which_gt0(delay_time);
      l[whichd] -= lS(dist, delay_time[whichd], eta[whichd], aux[whichd], aux2[whichd]);
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

  // Study IDs for independent shape parameters
  int<lower=1> study[narm_ipd + narm_agd_arm];
}
transformed data {
  // Exponential model indicator, 0 = exponential
  int<lower=0, upper=1> nonexp = (dist == 1 || dist == 4) ? 0 : 1;
  // Generalised Gamma model indicator, 1 = gengamma
  int<lower=0, upper=1> gengamma = dist >= 9 ? 1 : 0;

  // // Study ID for IPD studies
  // int<lower=1> ipd_study[ni_ipd] = study[1:ni_ipd];
  // // Study ID vector for AgD studies with expanded AgD integration points
  // int<lower=1> agd_arm_study[ni_ipd + nint * ni_agd_arm];

#include /include/transformed_data_common.stan

  // for (i in 1:ni_agd_arm) {
  //   agd_arm_study[((i-1)*nint + 1):(i*nint)] = rep_array(study[ni_ipd + i], nint);
  // }
}
parameters {
#include /include/parameters_common.stan

  // Auxiliary parameters for parametric model (typically shape)
  // Exponential model has shape = 1 so parameter is removed (zero dim)
  vector<lower=0>[(ns_ipd + ns_agd_arm)*nonexp] aux;
  // Second auxiliary parameter for generalised gamma
  vector<lower=0>[(ns_ipd + ns_agd_arm)*gengamma] aux2;
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
                     nonexp ? aux[study[ipd_arm]] : rep_vector(0, ni_ipd),
                     gengamma ? aux2[study[ipd_arm]] : rep_vector(0, ni_ipd));

  // -- AgD model (arm-based) --
  if (ni_agd_arm) {
    vector[nint * ni_agd_arm] eta_agd_arm_noRE = has_offset ?
              X_agd_arm * beta_tilde + offset_agd_arm :
              X_agd_arm * beta_tilde;

    if (nint > 1) { // -- If integration points are used --

      // Integration points are local variable only for efficiency
      vector[nint] eta_agd_arm_ii;
      vector[nint] log_L_ii;

      for (i in 1:ni_agd_arm) {
        eta_agd_arm_ii = eta_agd_arm_noRE[(1 + (i-1)*nint):(i*nint)];

        // Random effects
        if (RE && which_RE[narm_ipd + agd_arm_arm[i]]) eta_agd_arm_ii += f_delta[which_RE[narm_ipd + agd_arm_arm[i]]];

        // Average likelihood over integration points
        log_L_ii = loglik(dist,
                          agd_arm_time[i],
                          agd_arm_start_time[i],
                          agd_arm_delay_time[i],
                          agd_arm_status[i],
                          eta_agd_arm_ii,
                          nonexp ? aux[study[narm_ipd + agd_arm_arm[i]]] : 0,
                          gengamma ? aux2[study[narm_ipd + agd_arm_arm[i]]] : 0);

        log_L_agd_arm[i] = log_sum_exp(log_L_ii) - log(nint);
      }
    } else { // -- If no integration --
      for (i in 1:ni_agd_arm) {
        real eta_agd_arm = eta_agd_arm_noRE[i];

        // Random effects
        if (RE && which_RE[narm_ipd + agd_arm_arm[i]]) eta_agd_arm += f_delta[which_RE[narm_ipd + agd_arm_arm[i]]];

        log_L_agd_arm[i] = loglik(dist,
                                  agd_arm_time[i],
                                  agd_arm_start_time[i],
                                  agd_arm_delay_time[i],
                                  agd_arm_status[i],
                                  eta_agd_arm,
                                  nonexp ? aux[study[narm_ipd + agd_arm_arm[i]]] : 0,
                                  gengamma ? aux2[study[narm_ipd + agd_arm_arm[i]]] : 0);
      }
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

  vector[(dist != 1 && dist != 4 && dist != 6 && dist != 9) ? ns_ipd + ns_agd_arm : 0] shape;
  vector[dist == 6 ? ns_ipd + ns_agd_arm : 0] sdlog;  // lognormal sdlog
  vector[dist == 9 ? ns_ipd + ns_agd_arm : 0] sigma; // gengamma sigma
  vector[dist == 9 ? ns_ipd + ns_agd_arm : 0] k;  // gengamma k

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
