functions {
#include /include/prior_select.stan
#include /include/count_nonzero.stan

  //-- (log) Survival functions --
  real S(int dist, real y, real rate, real shape, int log_) {
    real S;

    if (dist == 1) { // Exponential
      if (log_ == 0) S = exp(-y * rate);
      else S = -y * rate;

    } else if (dist == 2) { // Weibull
      if (log_ == 0) {
        S = exp(-pow(y * rate, shape));
      } else {
        S = -pow(y * rate, shape);
      }

    } else if (dist == 3) { // Gompertz
      if (log_ == 0) {
        S = exp(-rate/shape * expm1(shape * y));
      } else {
        S = -rate/shape * expm1(shape * y);
      }
    }

    return S;
  }

  // -- (log) Hazard functions --
  real h(int dist, real y, real rate, real shape, int log_) {
    real h;

    if (dist == 1) { // Exponential
      if (log_ == 0) h = rate;
      else h = log(rate);

    } else if (dist == 2) { // Weibull
      if (log_ == 0) {
        h = shape * rate * pow(y * rate, shape - 1);
      } else {
        h = log(shape) + log(rate) + lmultiply(y * rate, shape - 1);
      }
    } else if (dist == 3) { // Gompertz
      if (log_ == 0) h = rate * exp(shape * y);
      else h = log(rate) + (shape * y);
    }

    return h;
  }

  // -- Likelihood with censoring and truncation --
  // vector lik(int dist, real time, real start_time, real delay_time, int status, vector rate, real shape) {
  //   vector[num_elements(rate)] l;
  //
  //   if (status == 0) { // Right censored
  //     l = S(dist, time, rate, shape, 0);
  //   } else if (status == 1) { // Observed
  //     l = S(dist, time, rate, shape, 0) .* h(dist, time, rate, shape, 0);
  //   } else if (status == 2) { // Left censored
  //     l = 1 - S(dist, time, rate, shape, 0);
  //   } else if (status == 3) { // Interval censored
  //     l = S(dist, start_time, rate, shape, 0) - S(dist, time, rate, shape, 0);
  //   }
  //
  //   // Left truncation
  //   if (delay_time > 0) {
  //     l = l ./ S(dist, delay_time, rate, shape, 0);
  //   }
  //
  //   return l;
  // }

  // -- Log likelihood with censoring and truncation --
  real loglik(int dist, real time, real start_time, real delay_time, int status, real rate, real shape) {
    real l;

    if (status == 0) { // Right censored
      l = S(dist, time, rate, shape, 1);
    } else if (status == 1) { // Observed
      l = S(dist, time, rate, shape, 1) + h(dist, time, rate, shape, 1);
    } else if (status == 2) { // Left censored
      l = log1m(S(dist, time, rate, shape, 0));
    } else if (status == 3) { // Interval censored
      l = log(S(dist, start_time, rate, shape, 0) - S(dist, time, rate, shape, 0));
    }

    // Left truncation
    if (delay_time > 0) {
      l -= S(dist, delay_time, rate, shape, 1);
    }

    return l;
  }
}
data {
#include /include/data_common.stan

  // Prior on shape parameter
  int<lower=0,upper=6> prior_aux_dist;
  real prior_aux_location;
  real<lower=0> prior_aux_scale;
  real<lower=0> prior_aux_df;

  // Select distribution
  int<lower=1, upper=3> dist;
  // 1 = Exponential
  // 2 = Weibull
  // 3 = Gompertz

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
  int<lower=0, upper=1> nonexp = dist == 1 ? 0 : 1;

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

  // Shape for parametric model
  // Exponential model has shape = 1 so parameter is removed (zero dim)
  vector<lower=0>[(ns_ipd + ns_agd_arm)*nonexp] shape;
}
transformed parameters {
  // Log likelihood contributions
  vector[ni_ipd] log_L_ipd;
  vector[ni_agd_arm] log_L_agd_arm;

#include /include/transformed_parameters_common.stan

  // -- IPD model --
  if (dist == 2) eta_ipd = eta_ipd ./ shape[study[ipd_arm]]; // Weibull linear predictor

  // Evaluate log likelihood
  for (i in 1:ni_ipd) {
    log_L_ipd[i] = loglik(dist,
                          ipd_time[i],
                          ipd_start_time[i],
                          ipd_delay_time[i],
                          ipd_status[i],
                          exp(eta_ipd[i]),
                          shape[study[ipd_arm[i]]]);
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

        // Weibull linear predictor
        if (dist == 2) eta_agd_arm_ii = eta_agd_arm_ii / shape[study[narm_ipd + agd_arm_arm[i]]];

        // Average likelihood over integration points
        // NOTE: Normalising term (dividing by nint) omitted as this is constant
        for (j in 1:nint) {
          log_L_ii[j] = loglik(dist,
                               agd_arm_time[i],
                               agd_arm_start_time[i],
                               agd_arm_delay_time[i],
                               agd_arm_status[i],
                               exp(eta_agd_arm_ii[j]),
                               shape[study[narm_ipd + agd_arm_arm[i]]]);
        }

        log_L_agd_arm[i] = log_sum_exp(log_L_ii);
      }
    } else { // -- If no integration --
      for (i in 1:ni_agd_arm) {
        real eta_agd_arm = eta_agd_arm_noRE[i];

        // Random effects
        if (RE && which_RE[narm_ipd + agd_arm_arm[i]]) eta_agd_arm += f_delta[which_RE[narm_ipd + agd_arm_arm[i]]];

        // Weibull linear predictor
        if (dist == 2) eta_agd_arm = eta_agd_arm / shape[study[narm_ipd + agd_arm_arm[i]]];

        log_L_agd_arm[i] = loglik(dist,
                                  agd_arm_time[i],
                                  agd_arm_start_time[i],
                                  agd_arm_delay_time[i],
                                  agd_arm_status[i],
                                  exp(eta_agd_arm),
                                  shape[study[narm_ipd + agd_arm_arm[i]]]);
      }
    }
  }
}
model {
#include /include/model_common.stan

  // -- Prior on shapes --
  if (nonexp) prior_select_lp(shape, prior_aux_dist, prior_aux_location, prior_aux_scale, prior_aux_df);

  // -- IPD likelihood --
  target += log_L_ipd;

  // -- AgD likelihood (arm-based) --
  target += log_L_agd_arm;
}
generated quantities {
  // Transform intercepts back to scales
  vector[ns_ipd + ns_agd_arm] scale = exp(mu);

#include /include/generated_quantities_theta_fitted.stan
#include /include/generated_quantities_common.stan

  // Log likelihood
  log_lik[1:ni_ipd] = log_L_ipd;
  log_lik[(ni_ipd + 1):(ni_ipd + ni_agd_arm)] = log_L_agd_arm - log(nint);

  // Residual deviance and fitted values

  // resdev[1:ni_ipd]
  // resdev[(ni_ipd + 1):(ni_ipd + ni_agd)]
  // fitted_ipd
  // fitted_agd_arm


}