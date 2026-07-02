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

// -- AgD model (regression coefficients) --
vector [ no_agd_regression               ? ni_agd_regression:0] eta_red;
vector [ no_agd_regression               ? ni_agd_regression:0] mu_red;
vector [(no_agd_regression && link != 1) ? ni_agd_regression:0] agd_regression_OVB_w; // non-canonical link function weights
array[no_agd_regression ? ns_agd_regression : 0] matrix[agd_regression_max_ncoef_inc, agd_regression_max_ncoef_omt] agd_regression_OVB_mat_omt; // (XI' XI)^{-1}XI' XO
array[no_agd_regression ? ns_agd_regression : 0] matrix[agd_regression_max_ncoef_inc, agd_regression_max_nrow     ] agd_regression_OVB_mat_hat; // (XI' XI)^{-1} XI'

if (no_agd_regression){
  int c_c = 0; // coef. counter
  int c_i = 0; // Included coef. counter
  int c_o = 0; // Omitted coef. counter
  int c_x = 0; // X_int rows counter
  for (i in 1:ns_agd_regression) {

    if(agd_regression_reduced_study[i]){

      // vector mdivide_left_spd( M,  N)  equals to inverse(M) * N
      // matrix crossprod(matrix x) equals to X'X

      agd_regression_OVB_mat_hat[i][1:agd_regression_ncoef_inc[i], 1:agd_regression_nx[i] ]  =
      mdivide_left_spd( crossprod(X_agd_regression_int[ (c_x+1):(c_x+agd_regression_nx[i]) ,XI_col_vec[(c_i+1):(c_i+agd_regression_ncoef_inc[i])] ]) ,
      (X_agd_regression_int[ (c_x+1):(c_x+agd_regression_nx[i]) ,XI_col_vec[(c_i+1):(c_i+agd_regression_ncoef_inc[i])] ])' );

      agd_regression_OVB_mat_omt[i][1:agd_regression_ncoef_inc[i], 1:agd_regression_ncoef_omt[i] ]  =
      agd_regression_OVB_mat_hat[i][1:agd_regression_ncoef_inc[i], 1:agd_regression_nx[i] ] *
      (X_agd_regression_int[ (c_x+1):(c_x+agd_regression_nx[i]) ,XO_col_vec[(c_o+1):(c_o+agd_regression_ncoef_omt[i])] ]);

      eta_red[(c_x+1):(c_x+agd_regression_nx[i])] =
      X_agd_regression_int[(c_x+1):(c_x+agd_regression_nx[i]) ,XI_col_vec[(c_i+1):(c_i+agd_regression_ncoef_inc[i])]] *
      agd_regression_est[(c_c+1):(c_c+agd_regression_ncoef[i])];

      if (link == 1){ // logit link
       mu_red[(c_x+1):(c_x+agd_regression_nx[i])]  =  inv_logit(eta_red[(c_x+1):(c_x+agd_regression_nx[i])]);
      }else if (link == 2){ // probit link
       mu_red[(c_x+1):(c_x+agd_regression_nx[i])]  =  Phi(eta_red[(c_x+1):(c_x+agd_regression_nx[i])]);
       agd_regression_OVB_w[(c_x+1):(c_x+agd_regression_nx[i])] =
       exp(std_normal_lpdf(eta_red[(c_x+1):(c_x+agd_regression_nx[i])])) ./
       (mu_red[(c_x+1):(c_x+agd_regression_nx[i])].*(1-mu_red[(c_x+1):(c_x+agd_regression_nx[i])]));
      }else if (link == 3){ // cloglog link
       mu_red[(c_x+1):(c_x+agd_regression_nx[i])]  =  inv_cloglog(eta_red[(c_x+1):(c_x+agd_regression_nx[i])]);
       agd_regression_OVB_w[(c_x+1):(c_x+agd_regression_nx[i])] =
       exp(mu_red[(c_x+1):(c_x+agd_regression_nx[i])]) ./
       mu_red[(c_x+1):(c_x+agd_regression_nx[i])];
      }
    }

    c_c += agd_regression_ncoef[i];
    c_i += agd_regression_ncoef_inc[i];
    c_o += agd_regression_ncoef_omt[i];
    c_x += agd_regression_nx[i];
  }
}


}
parameters {
#include /include/parameters_common.stan

}
transformed parameters {
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
    vector[nint_max * ni_agd_arm] eta_agd_arm_noRE = has_offset ?
      X_agd_arm * beta_tilde + offset_agd_arm :
      X_agd_arm * beta_tilde;

    // Baseline risk meta-regression
    if (brmr_n_col > 0) {
      // Subtracting 1 from the centred baseline risk here, as the associated
      // beta was already added once to the linear predictor by
      // `X_agd_arm * beta_tilde`
      eta_agd_arm_noRE += (X_agd_arm[,1:totns] * mu - xbar_mu - 1) .* (X_agd_arm[,brmr_col] * beta_tilde[brmr_col]);
    }

    if (nint_max > 1) { // -- If integration points are used --

      if (class_effects) {
        for (i in 1:ni_agd_arm) {
          if (agd_arm_trt[i] > 1 && which_CE[agd_arm_trt[i] - 1]) {
            eta_agd_arm_noRE[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)] += f_class[which_fclass[agd_arm_trt[i] - 1]];
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

        for (i in 1:ni_agd_arm) {
          theta_agd_arm_bar[i] = mean(theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint)]);
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

        for (i in 1:ni_agd_arm) {
          theta_agd_arm_bar[i] = mean(theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint)]);
        }
      }
    } else { // -- If no integration --
      if (RE) {

        // Add class effects contribution to the linear predictor
        if (class_effects) {
          for (i in 1:ni_agd_arm) {
            if (agd_arm_trt[i] > 1 && which_CE[agd_arm_trt[i] - 1]) {
              eta_agd_arm_noRE[i] += f_class[which_fclass[agd_arm_trt[i] - 1]];
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

        // Add class effects contribution to the linear predictor
        if (class_effects) {
          for (i in 1:ni_agd_arm) {
            if (agd_arm_trt[i] > 1 && which_CE[agd_arm_trt[i] - 1]) {
              eta_agd_arm_noRE[i] += f_class[which_fclass[agd_arm_trt[i] - 1]];
            }
          }
        }

        if (link == 1) // logit link
          theta_agd_arm_bar = inv_logit(eta_agd_arm_noRE);
        else if (link == 2) // probit link
          theta_agd_arm_bar = Phi(eta_agd_arm_noRE);
        else if (link == 3) // cloglog link
          theta_agd_arm_bar = inv_cloglog(eta_agd_arm_noRE);
      }
    }
  }

// -- AgD model (regression coefficients) --
vector [no_agd_regression ? ni_agd_regression:0] mu_ful;
vector [no_agd_regression ? ni_agd_regression:0] lp_err; // linear predictor mismatch
vector [no_agd_regression ? ni_agd_regression:0] mu_err; // mu mismatch

if (nc_agd_regression) {

  if (sum(agd_regression_reduced_study)){

    int c_c = 0; // coef. counter
    int c_i = 0; // Included coef. counter
    int c_o = 0; // Omitted coef. counter
    int c_x = 0; // X_int rows counter

    for (i in 1:ns_agd_regression) {
      // OVB adjustment
      if (agd_regression_reduced_study[i]){

        if (link == 1){ // logit link

          mu_ful[(c_x+1):(c_x+agd_regression_nx[i])] =  inv_logit(X_agd_regression_int[ (c_x+1):(c_x+agd_regression_nx[i]),] * allbeta) ;

          mu_err[(c_x+1):(c_x+agd_regression_nx[i])] =
          inv_logit(mu_red[ (c_x+1):(c_x+agd_regression_nx[i])]) -
          inv_logit(mu_ful[ (c_x+1):(c_x+agd_regression_nx[i])]);

          lp_err[(c_x+1):(c_x+agd_regression_nx[i])] =
          eta_red[(c_x+1):(c_x+agd_regression_nx[i])] -
          logit(mu_red[(c_x+1):(c_x+agd_regression_nx[i])] - mu_err[(c_x+1):(c_x+agd_regression_nx[i])]);

          eta_agd_regression[(c_c+1):(c_c+agd_regression_ncoef[i])] =
            allbeta[XI_col_vec[(c_i+1):(c_i+agd_regression_ncoef_inc[i])]] +
            block(agd_regression_OVB_mat_omt[i], 1, 1,agd_regression_ncoef_inc[i] ,agd_regression_ncoef_omt[i] ) *
            allbeta[XO_col_vec[(c_o+1):(c_o+agd_regression_ncoef_omt[i])]] +
            block(agd_regression_OVB_mat_hat[i]  , 1, 1,agd_regression_ncoef_inc[i] ,agd_regression_nx[i] ) *
            lp_err[(c_x+1):(c_x+agd_regression_nx[i])];


        }else if (link == 2){ // probit link

          mu_ful[(c_x+1):(c_x+agd_regression_nx[i])] =  Phi(X_agd_regression_int[ (c_x+1):(c_x+agd_regression_nx[i]),] * allbeta) ;

          mu_err[(c_x+1):(c_x+agd_regression_nx[i])] =
          Phi(mu_red[ (c_x+1):(c_x+agd_regression_nx[i])]) -
          Phi(mu_ful[ (c_x+1):(c_x+agd_regression_nx[i])]);

          lp_err[(c_x+1):(c_x+agd_regression_nx[i])] =
          eta_red[(c_x+1):(c_x+agd_regression_nx[i])] -
          inv_Phi(mu_red[(c_x+1):(c_x+agd_regression_nx[i])] - mu_err[(c_x+1):(c_x+agd_regression_nx[i])]);

          eta_agd_regression[(c_c+1):(c_c+agd_regression_ncoef[i])] =
            allbeta[XI_col_vec[(c_i+1):(c_i+agd_regression_ncoef_inc[i])]] +
            block(agd_regression_OVB_mat_omt[i], 1, 1,agd_regression_ncoef_inc[i] ,agd_regression_ncoef_omt[i] ) *
            allbeta[XO_col_vec[(c_o+1):(c_o+agd_regression_ncoef_omt[i])]] +
            block(agd_regression_OVB_mat_hat[i]  , 1, 1,agd_regression_ncoef_inc[i] ,agd_regression_nx[i] ) *
            lp_err[(c_x+1):(c_x+agd_regression_nx[i])];

        }else if (link == 3){ // cloglog link

          mu_ful[(c_x+1):(c_x+agd_regression_nx[i])] =  inv_cloglog(X_agd_regression_int[ (c_x+1):(c_x+agd_regression_nx[i]),] * allbeta) ;

          mu_err[(c_x+1):(c_x+agd_regression_nx[i])] =
          inv_cloglog(mu_red[ (c_x+1):(c_x+agd_regression_nx[i])]) -
          inv_cloglog(mu_ful[ (c_x+1):(c_x+agd_regression_nx[i])]);

          lp_err[(c_x+1):(c_x+agd_regression_nx[i])] =
          eta_red[(c_x+1):(c_x+agd_regression_nx[i])] -
          log(-log1m(mu_red[(c_x+1):(c_x+agd_regression_nx[i])] - mu_err[(c_x+1):(c_x+agd_regression_nx[i])]));

          eta_agd_regression[(c_c+1):(c_c+agd_regression_ncoef[i])] =
            allbeta[XI_col_vec[(c_i+1):(c_i+agd_regression_ncoef_inc[i])]] +
            block(agd_regression_OVB_mat_omt[i], 1, 1,agd_regression_ncoef_inc[i] ,agd_regression_ncoef_omt[i] ) *
            allbeta[XO_col_vec[(c_o+1):(c_o+agd_regression_ncoef_omt[i])]] +
            block(agd_regression_OVB_mat_hat[i]  , 1, 1,agd_regression_ncoef_inc[i] ,agd_regression_nx[i] ) *
            lp_err[(c_x+1):(c_x+agd_regression_nx[i])];

        }
      }else{
        eta_agd_regression[ (c_c+1):(c_c+agd_regression_ncoef[i]) ] = X_agd_regression_no_QR[ (c_c+1):(c_c+agd_regression_ncoef[i]), ] * allbeta;
      }
      c_c += agd_regression_ncoef[i];
      c_i += agd_regression_ncoef_inc[i];
      c_o += agd_regression_ncoef_omt[i];
      c_x += agd_regression_nx[i];
    }
  }else{
    eta_agd_regression = X_agd_regression * beta_tilde;
  }

  if (RE) {
    for (i in 1:nc_agd_regression) {
      if (which_RE[narm_ipd + narm_agd_arm + ni_agd_contrast + i])
        eta_agd_regression[i] = eta_agd_regression[i] + f_delta[which_RE[narm_ipd + narm_agd_arm + ni_agd_contrast + i]];
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
  agd_arm_r ~ binomial(agd_arm_n, theta_agd_arm_bar);

  // -- AgD regression soft constraints --
  if(no_agd_regression){
    int c_x = 0;
    int c_i = 0;
    for (i in 1:ns_agd_regression) {
      if(agd_regression_reduced_study[i]){

          if (link == 1){ // logit link (canonical link function)

          (
            ((X_agd_regression_int[(c_x+1):(c_x+agd_regression_nx[i]) , XI_col_vec[(c_i+1):(c_i+agd_regression_ncoef_inc[i])]])' *
            mu_err[(c_x+1):(c_x+agd_regression_nx[i])]) / agd_regression_nx[i]
            ) ~ normal( 0 , 0.01) ;

          } else { // (non-canonical link function)

            (
              ((X_agd_regression_int[(c_x+1):(c_x+agd_regression_nx[i]) , XI_col_vec[(c_i+1):(c_i+agd_regression_ncoef_inc[i])]])' *
              (mu_err[(c_x+1):(c_x+agd_regression_nx[i])] .* agd_regression_OVB_w)) / agd_regression_nx[i]
              ) ~ normal( 0 , 0.01) ;

          }

      }
      c_i += agd_regression_ncoef_inc[i];
      c_x += agd_regression_nx[i];
    }
  }

}
generated quantities {
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
    log_lik[ni_ipd + i] = binomial_lpmf(agd_arm_r[i] | agd_arm_n[i], theta_agd_arm_bar[i]);
    resdev[ni_ipd + i] = 2 *
      ((agd_arm_r[i] > 0 ?
         lmultiply(agd_arm_r[i],
                   agd_arm_r[i] / (agd_arm_n[i] * theta_agd_arm_bar[i])) : 0) +
       (agd_arm_r[i] < agd_arm_n[i] ?
         lmultiply(agd_arm_n[i] - agd_arm_r[i],
                  (agd_arm_n[i] - agd_arm_r[i]) / (agd_arm_n[i] - agd_arm_n[i] * theta_agd_arm_bar[i])) : 0));
    fitted_agd_arm[i] = agd_arm_n[i] * theta_agd_arm_bar[i];
  }

}
