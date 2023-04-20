functions {
#include /include/prior_select.stan
#include /include/count_nonzero.stan

  // -- Generalised Gamma ldpf --
  real gengamma_lpdf(real y, real mu, real sigma, real k) {
    // Parameterisation of Lawless
	  real Q = pow(k, -0.5);
	  real w = Q * (log(y) - mu) / sigma;
	  return -log(sigma) - log(y) - 0.5 * log(k) * (1 - 2 * k) + k * (w - exp(w)) - lgamma(k);
  }

  //-- log Survival functions --
  real lS(int dist, real y, real eta, real aux, real aux2) {
    // aux is shape, except for lognormal sdlog, gengamma sigma
    // aux2 is only for gengamma k

    real lS;

    if (dist == 1) { // Exponential
      lS = -y * exp(eta);
    } else if (dist == 2) { // Weibull
      lS = -pow(y, aux) * exp(eta);
    } else if (dist == 3) { // Gompertz
      lS = -exp(eta)/aux * expm1(aux * y);
    } else if (dist == 4) { // Exponential AFT
      lS = -y * exp(-eta);
    } else if (dist == 5) { // Weibull AFT
      lS = -pow(y, aux) * exp(-aux * eta);
    } else if (dist == 6) { // log Normal
      // aux is sdlog
      lS = log1m(Phi((log(y) - eta) / aux));
    } else if (dist == 7) { // log logistic
      lS = -log1p(pow(y / exp(eta), aux));
    } else if (dist == 8) { // Gamma
      lS = gamma_lccdf(y | aux, exp(-eta));
    } else if (dist == 9) { // Generalised Gamma
      real Q = pow(aux2, -0.5);
      real w = exp(Q * (log(y) - eta) / aux) * aux2;
      lS = log1m(gamma_p(aux2, w));
    }

    return lS;
  }

  // -- log Hazard functions --
  real lh(int dist, real y, real eta, real aux, real aux2) {
    // aux is shape, except for lognormal sdlog, gengamma scale
    // aux2 is only for gengamma shape

    real lh;

    if (dist == 1) { // Exponential
      lh = eta;
    } else if (dist == 2) { // Weibull
      lh = log(aux) + eta + lmultiply(y, aux - 1);
    } else if (dist == 3) { // Gompertz
      lh = eta + (aux * y);
    } else if (dist == 4) { // Exponential AFT
      lh = -eta;
    } else if (dist == 5) { // Weibull AFT
      lh = log(aux) - (aux * eta) + lmultiply(y, aux - 1);
    } else if (dist == 6) { // log Normal
      // aux is sdlog
      lh = lognormal_lpdf(y | eta, aux) - log1m(Phi((log(y) - eta) / aux));
    } else if (dist == 7) { // log logistic
      lh = log(aux) - eta + (aux - 1)*(log(y) - eta) - log1p(pow(y / exp(eta), aux));
    } else if (dist == 8) { // Gamma
      lh = gamma_lpdf(y | aux, exp(-eta)) - gamma_lccdf(y | aux, exp(-eta));
    } else if (dist == 9) { // Generalised Gamma
      // Not used, lpdf used directly
    }

    return lh;
  }

  // -- Log likelihood with censoring and truncation --
  real loglik(int dist, real time, real start_time, real delay_time, int status, real eta, real aux, real aux2) {
    real l;

    if (status == 0) { // Right censored
      l = lS(dist, time, eta, aux, aux2);
    } else if (status == 1) { // Observed
      // Make (generalised) Gamma models more efficient by using ldpf directly
      if (dist == 8) { // Gamma
        l = gamma_lpdf(time | aux, exp(-eta));
      } if (dist == 9) { // Gen Gamma
        l = gengamma_lpdf(time | eta, aux, aux2);
      } else {
        l = lS(dist, time, eta, aux, aux2) + lh(dist, time, eta, aux, aux2);
      }
    } else if (status == 2) { // Left censored
      l = log1m_exp(lS(dist, time, eta, aux, aux2));
    } else if (status == 3) { // Interval censored
      l = log_diff_exp(lS(dist, start_time, eta, aux, aux2), lS(dist, time, eta, aux, aux2));
    }

    // Left truncation
    if (delay_time > 0) {
      l -= lS(dist, delay_time, eta, aux, aux2);
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
  real<lower=0> ipd_time[ni_ipd];
  real<lower=0> ipd_start_time[ni_ipd];
  real<lower=0> ipd_delay_time[ni_ipd];
  int<lower=0, upper=3> ipd_status[ni_ipd];

  real<lower=0> agd_arm_time[ni_agd_arm];
  real<lower=0> agd_arm_start_time[ni_agd_arm];
  real<lower=0> agd_arm_delay_time[ni_agd_arm];
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
  for (i in 1:ni_ipd) {
    log_L_ipd[i] = loglik(dist,
                          ipd_time[i],
                          ipd_start_time[i],
                          ipd_delay_time[i],
                          ipd_status[i],
                          eta_ipd[i],
                          nonexp ? aux[study[ipd_arm[i]]] : 0,
                          gengamma ? aux2[study[ipd_arm[i]]] : 0);
  }

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
        // NOTE: Normalising term (dividing by nint) omitted as this is constant
        for (j in 1:nint) {
          log_L_ii[j] = loglik(dist,
                               agd_arm_time[i],
                               agd_arm_start_time[i],
                               agd_arm_delay_time[i],
                               agd_arm_status[i],
                               eta_agd_arm_ii[j],
                               nonexp ? aux[study[narm_ipd + agd_arm_arm[i]]] : 0,
                               gengamma ? aux2[study[narm_ipd + agd_arm_arm[i]]] : 0);
        }

        log_L_agd_arm[i] = log_sum_exp(log_L_ii);
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
  log_lik[1:ni_ipd] = log_L_ipd;
  log_lik[(ni_ipd + 1):(ni_ipd + ni_agd_arm)] = log_L_agd_arm - log(nint);

  // Residual deviance
  // NOTE: This is actually just the deviance, which is fine for relative fit model comparison
  resdev[1:(ni_ipd + ni_agd_arm)] = -2 * log_lik[1:(ni_ipd + ni_agd_arm)];

  // Fitted values not implemented

  // Rename parameters
  if (dist != 1 && dist != 4 && dist != 6 && dist != 9) shape = aux;
  if (dist == 6) sdlog = aux;
  if (dist == 9) {
    sigma = aux;
    k = aux2;
  }

}
