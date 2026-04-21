functions {
#include /include/prior_select.stan
#include /include/count_nonzero.stan
}
data {
#include /include/data_common.stan

  // -- AgD regression coefficients --
  int<lower=0> agd_regression_glm_method;
  vector[no_agd_regression ? nc_agd_regression : 0] agd_regression_OVB_GLM; //  GLM OVB adjustment

  // Outcomes
  array[ni_ipd] int<lower=0, upper=1> ipd_r;
  array[ni_agd_arm] int<lower=0> agd_arm_n;
  array[ni_agd_arm] int<lower=0> agd_arm_r;
}
transformed data {
#include /include/transformed_data_common.stan
}
parameters {
#include /include/parameters_common.stan
real<lower=0> err_sd;
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
vector[ni_agd_regression] err;
matrix[ni_agd_regression,nX] err_mat;
if (nc_agd_regression) {

  if (sum(agd_regression_reduced_study)){
    // vector [nX] allbeta_OVB;
    int c_c = 0; // coef. counter
    int c_i = 0; // Included coef. counter
    // int c_o = 0; // Omitted coef. counter
    int c_x = 0; // X_int rows counter
    for (i in 1:ns_agd_regression) {
      // allbeta_OVB = allbeta;
      // OVB adjustment
      if(agd_regression_reduced_study[i]){
        if (link == 1){ // logit link
          // eta_agd_regression[ (c_c+1):(c_c+agd_regression_ncoef[i]) ] = logit(mean(inv_logit(       X_agd_regression_int[(c_x+1):(c_x+agd_regression_nx[i]),]*allbeta)))  - agd_regression_OVB_GLM[(c_c+1):(c_c+agd_regression_ncoef[i])];

          # Calculate differences between the predicted responses from full and reduced models
          err[(c_x+1):(c_x+agd_regression_nx[i])] =
          inv_logit(X_agd_regression_int[ (c_x+1):(c_x+agd_regression_nx[i]) ,XI_col_vec[(c_i+1):(c_i+agd_regression_ncoef_inc[i])] ] * agd_regression_est[ (c_c+1):(c_c+agd_regression_ncoef[i]) ]) -
          inv_logit(X_agd_regression_int[ (c_x+1):(c_x+agd_regression_nx[i]) ,                                                      ] * allbeta) ;

          # Calculate differences between predicted responses between full and reduced models, multiplied by each included covariate
          err_mat[              (c_x+1):(c_x+agd_regression_nx[i]), XI_col_vec[(c_i+1):(c_i+agd_regression_ncoef_inc[i])] ] =
          X_agd_regression_int[ (c_x+1):(c_x+agd_regression_nx[i]) ,XI_col_vec[(c_i+1):(c_i+agd_regression_ncoef_inc[i])] ] .*
                 rep_matrix( err[(c_x+1):(c_x+agd_regression_nx[i])] ,  agd_regression_ncoef_inc[i] ) ;

          //  B = inverse(D'D)*D'X_O:
          //    - If M = D'D  and N = D'X_O
          //    - vector mdivide_left_spd( M,  N)  equals to inverse(M) * N
          //    - also, matrix crossprod(matrix x) equals to X'X
          eta_agd_regression[ (c_c+1):(c_c+agd_regression_ncoef[i]) ] =
          mdivide_left_spd( crossprod(X_agd_regression_int[ (c_x+1):(c_x+agd_regression_nx[i]) ,XI_col_vec[(c_i+1):(c_i+agd_regression_ncoef_inc[i])] ]) ,
          (X_agd_regression_int[ (c_x+1):(c_x+agd_regression_nx[i]) ,XI_col_vec[(c_i+1):(c_i+agd_regression_ncoef_inc[i])] ])' ) *
          logit(inv_logit( X_agd_regression_int[(c_x+1):(c_x+agd_regression_nx[i]),]*allbeta)  +  err[(c_x+1):(c_x+agd_regression_nx[i])]) ;

        }else if (link == 2){ // probit link
          eta_agd_regression[ (c_c+1):(c_c+agd_regression_ncoef[i]) ] = inv_Phi(mean(Phi(           X_agd_regression_int[(c_x+1):(c_x+agd_regression_nx[i]),]*allbeta)))  - agd_regression_OVB_GLM[(c_c+1):(c_c+agd_regression_ncoef[i])];
        }else if (link == 3){ // cloglog link
          eta_agd_regression[ (c_c+1):(c_c+agd_regression_ncoef[i]) ] = log(-log1m(mean(inv_cloglog(X_agd_regression_int[(c_x+1):(c_x+agd_regression_nx[i]),]*allbeta)))) - agd_regression_OVB_GLM[(c_c+1):(c_c+agd_regression_ncoef[i])];
        }
      }else{
        eta_agd_regression[ (c_c+1):(c_c+agd_regression_ncoef[i]) ] = X_agd_regression_no_QR[ (c_c+1):(c_c+agd_regression_ncoef[i]), ] * allbeta;
      }
      c_c += agd_regression_ncoef[i];
      c_i += agd_regression_ncoef_inc[i];
      // c_o += agd_regression_ncoef_omt[i];
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

  if(nc_agd_regression){
    int c_x = 0;
    int c_i = 0;
    for (i in 1:ns_agd_regression) {
      if(agd_regression_reduced_study[i]){
        if(agd_regression_glm_method == 1 ){
          mean(err[(c_x+1):(c_x+agd_regression_nx[i])]) ~ normal( 0 , err_sd);
        } else if(agd_regression_glm_method == 2 ){
              (err[(c_x+1):(c_x+agd_regression_nx[i])]) ~ normal( 0 , err_sd);
        } else if(agd_regression_glm_method == 3 ){
          for (j in 1:agd_regression_ncoef_inc[i]) {
            mean((err_mat[(c_x+1):(c_x+agd_regression_nx[i]), XI_col_vec[(c_i+1):(c_i+agd_regression_ncoef_inc[i])] ])[,j])~ normal( 0 , err_sd);
          }
        } else if(agd_regression_glm_method == 4 ){
          for (j in 1:agd_regression_ncoef_inc[i]) {
                ((err_mat[(c_x+1):(c_x+agd_regression_nx[i]), XI_col_vec[(c_i+1):(c_i+agd_regression_ncoef_inc[i])] ])[,j])~ normal( 0 , err_sd);
          }
        }
      }
      c_i += agd_regression_ncoef_inc[i];
      c_x += agd_regression_nx[i];
    }
    err_sd ~ cauchy(0, 2.5);
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
