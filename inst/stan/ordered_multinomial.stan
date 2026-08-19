functions {
#include /include/prior_select.stan
#include /include/count_nonzero.stan
}
data {
#include /include/data_common.stan

  // Outcomes
  int<lower=2> ncat;

  array[ni_ipd] int<lower=1, upper=ncat> ipd_r;
  array[ni_ipd, ncat] int<lower=0, upper=ncat> ipd_cat;  // Category specs, left-aligned, padded with zeros
  array[ni_ipd] int<lower=2, upper=ncat> ipd_ncat;  // Number of categories observed

  array[ni_agd_arm, ncat] int<lower=0> agd_arm_r;
  vector[ni_agd_arm] agd_arm_n; // AgD arm sample sizes
  array[ni_agd_arm, ncat] int<lower=0, upper=ncat> agd_arm_cat;  // Category specs, left-aligned, padded with zeros
  array[ni_agd_arm] int<lower=2, upper=ncat> agd_arm_ncat;  // Number of categories observed

  // -- AgD regression coefficients --
  matrix[ns_agd_regression ? nc_agd_regression : 0, ns_agd_regression ? (ncat-1) : 0] X_agd_regression_cc; // Complementary design matrix for cc
  array [no_agd_regression ? nc_agd_regression : 0] int agd_regression_est_map_cc; // Map each position in the coefficient vector with intercepts to its position in the vector without intercepts (except study baseline)

  // Prior on differences between cutpoints
  int<lower=0,upper=6> prior_aux_dist;
  real prior_aux_location;
  real<lower=0> prior_aux_scale;
  real<lower=0> prior_aux_df;

}
transformed data {
  array[ni_ipd] vector[ncat] theta_ipd0 = rep_array(rep_vector(0, ncat), ni_ipd);
  array[ni_agd_arm] vector[ncat] theta_agd_arm_bar0 = rep_array(rep_vector(0, ncat), ni_agd_arm);
  // matrix[ni_agd_arm * n_int_thin, ncat] theta_bar_cum_agd_arm0 = rep_matrix(0, ni_agd_arm * n_int_thin, ncat);
  int ncat_mid = (ncat - 2) / 2 + 1; // middle cat
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

      int tmp_idx [agd_regression_max_ncoef];
      int tmp_idx_n = 0;
      for (j in 1:agd_regression_ncoef[i] ) {
        // index rows that are NOT related to cutpoints, except the smallest one (study baseline)
        if (sum(row(X_agd_regression_cc[(c_c+1):(c_c+agd_regression_ncoef[i]), 2:(ncat-1)], j)) == 0) {
          tmp_idx_n += 1;
          tmp_idx[tmp_idx_n] = j;
        }
      }

      eta_red[(c_x+1):(c_x+agd_regression_nx[i])] =
      X_agd_regression_int[(c_x+1):(c_x+agd_regression_nx[i]) ,XI_col_vec[(c_i+1):(c_i+agd_regression_ncoef_inc[i])]] *
      agd_regression_est[(c_c+1):(c_c+agd_regression_ncoef[i])][tmp_idx[1:tmp_idx_n]];
      // agd_regression_est[(c_c+1):(c_c+agd_regression_ncoef[i])];

      if (link == 1){ // logit link
       mu_red[(c_x+1):(c_x+agd_regression_nx[i])]  =  1 - inv_logit(eta_red[(c_x+1):(c_x+agd_regression_nx[i])]);
      }else if (link == 2){ // probit link
       mu_red[(c_x+1):(c_x+agd_regression_nx[i])]  =  1 - Phi(eta_red[(c_x+1):(c_x+agd_regression_nx[i])]);
       agd_regression_OVB_w[(c_x+1):(c_x+agd_regression_nx[i])] =
       exp(std_normal_lpdf(eta_red[(c_x+1):(c_x+agd_regression_nx[i])])) ./
       (mu_red[(c_x+1):(c_x+agd_regression_nx[i])].*(1-mu_red[(c_x+1):(c_x+agd_regression_nx[i])]));
      }else if (link == 3){ // cloglog link
       mu_red[(c_x+1):(c_x+agd_regression_nx[i])]  =  1 - inv_cloglog(eta_red[(c_x+1):(c_x+agd_regression_nx[i])]);
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

  // Ordered cutoffs on underlying probit-PASI scale
  // "Fixed effect" cutoffs, the same across trials
  positive_ordered[ncat - 2] f_cc;
}
transformed parameters {
  vector[ncat - 1] cc;

  array[ni_ipd] vector[ncat] theta_ipd; // IPD transformed predictor

  matrix[nint_max > 1 ? nint * ni_agd_arm : 0, nint_max > 1 ? ncat - 1 : 0] theta_agd_arm_ii; // Use these as q_ii intermediates
  array[ni_agd_arm] vector[ncat] q_agd_arm_bar; // AgD arm transformed predictor
  array[ni_agd_arm] vector[ncat] theta_agd_arm_bar; // AgD arm transformed predictor

#include /include/transformed_parameters_common.stan

  cc[1] = 0;
  if (ncat > 2) cc[2:] = f_cc;

  // Set predictors to zero for all missing categories (will drop out of likelihood)
  theta_ipd = theta_ipd0;
  theta_agd_arm_bar = theta_agd_arm_bar0;

  // -- IPD model --
  // Is this only necessary if link > 2? Since ordered_(logistic|probit) are available
  for (i in 1:ni_ipd) {
    vector[ipd_ncat[i] - 1] q_temp;
    for (k in 1:(ipd_ncat[i] - 1)) {
      if (link == 1) // logit link
        q_temp[k] = inv_logit(eta_ipd[i] - cc[ipd_cat[i, k+1]-1]);
      else if (link == 2) // probit link
        q_temp[k] = Phi(eta_ipd[i] - cc[ipd_cat[i, k+1]-1]);
      else if (link == 3) // cloglog link
        q_temp[k] = inv_cloglog(eta_ipd[i] - cc[ipd_cat[i, k+1]-1]);
    }

    // Category 1
    theta_ipd[i, 1] = 1 - q_temp[1];

    // Categories 2:(ipd_ncat - 1)
    for (k in 2:(ipd_ncat[i] - 1)) {
      // Store predictor in actual category column, rather than left-aligned
      theta_ipd[i, ipd_cat[i, k]] = q_temp[k - 1] - q_temp[k];
    }

    // Category ipd_ncat
    theta_ipd[i, ipd_cat[i, ipd_ncat[i]]] = q_temp[ipd_ncat[i] - 1];
  }

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

    // Add class effects contribution to the linear predictor
    if (class_effects) {
      for (i in 1:ni_agd_arm) {
        if (agd_arm_trt[i] > 1 && which_CE[agd_arm_trt[i] - 1]) {
          eta_agd_arm_noRE[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)] += f_class[which_fclass[agd_arm_trt[i] - 1]];
        }
      }
    }

    if (nint_max > 1) { // -- If integration points are used --

      if (RE) {

        vector[nint] eta_agd_arm_RE;

        if (link == 1) { // logit link
          for (i in 1:ni_agd_arm) {
            if (which_RE[narm_ipd + i])
              eta_agd_arm_RE = eta_agd_arm_noRE[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)] + f_delta[which_RE[narm_ipd + i]];
            else
              eta_agd_arm_RE = eta_agd_arm_noRE[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)];

            for (k in 1:(agd_arm_ncat[i] - 1)) {
              theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint), k] = inv_logit(eta_agd_arm_RE - cc[agd_arm_cat[i, k+1]-1]);
            }
          }
        } else if (link == 2) { // probit link
          for (i in 1:ni_agd_arm) {
            if (which_RE[narm_ipd + i])
              eta_agd_arm_RE = eta_agd_arm_noRE[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)] + f_delta[which_RE[narm_ipd + i]];
            else
              eta_agd_arm_RE = eta_agd_arm_noRE[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)];

            for (k in 1:(agd_arm_ncat[i] - 1)) {
              theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint), k] = Phi(eta_agd_arm_RE - cc[agd_arm_cat[i, k+1]-1]);
            }
          }
        } else if (link == 3) { // cloglog link
          for (i in 1:ni_agd_arm) {
            if (which_RE[narm_ipd + i])
              eta_agd_arm_RE = eta_agd_arm_noRE[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)] + f_delta[which_RE[narm_ipd + i]];
            else
              eta_agd_arm_RE = eta_agd_arm_noRE[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)];

            for (k in 1:(agd_arm_ncat[i] - 1)) {
              theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint), k] = inv_cloglog(eta_agd_arm_RE - cc[agd_arm_cat[i, k+1]-1]);
            }
          }
        }

      } else {

        if (link == 1) { // logit link
          for (i in 1:ni_agd_arm) {
            for (k in 1:(agd_arm_ncat[i] - 1)) {
              theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint), k] = inv_logit(eta_agd_arm_noRE[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)] - cc[agd_arm_cat[i, k+1]-1]);
            }
          }
        } else if (link == 2) { // probit link
          for (i in 1:ni_agd_arm) {
            for (k in 1:(agd_arm_ncat[i] - 1)) {
              theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint), k] = Phi(eta_agd_arm_noRE[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)] - cc[agd_arm_cat[i, k+1]-1]);
            }
          }
        } else if (link == 3) { // cloglog link
          for (i in 1:ni_agd_arm) {
            for (k in 1:(agd_arm_ncat[i] - 1)) {
              theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint), k] = inv_cloglog(eta_agd_arm_noRE[(1 + (i-1)*nint_max):((i-1)*nint_max + nint)] - cc[agd_arm_cat[i, k+1]-1]);
            }
          }
        }

      }

      for (i in 1:ni_agd_arm) {
        for (k in 1:(agd_arm_ncat[i] - 1)) {
          q_agd_arm_bar[i, k] = mean(theta_agd_arm_ii[(1 + (i-1)*nint):(i*nint), k]);
        }
      }

    } else { // -- If no integration --

      if (RE) {

        real eta_agd_arm_RE;

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
              eta_agd_arm_RE = eta_agd_arm_noRE[i] + f_delta[which_RE[narm_ipd + i]];
            else
              eta_agd_arm_RE = eta_agd_arm_noRE[i];

            for (k in 1:(agd_arm_ncat[i] - 1)) {
              q_agd_arm_bar[i, k] = inv_logit(eta_agd_arm_RE - cc[agd_arm_cat[i, k+1]-1]);
            }
          }
        } else if (link == 2) { // probit link
          for (i in 1:ni_agd_arm) {
            if (which_RE[narm_ipd + i])
              eta_agd_arm_RE = eta_agd_arm_noRE[i] + f_delta[which_RE[narm_ipd + i]];
            else
              eta_agd_arm_RE = eta_agd_arm_noRE[i];

            for (k in 1:(agd_arm_ncat[i] - 1)) {
              q_agd_arm_bar[i, k] = Phi(eta_agd_arm_RE - cc[agd_arm_cat[i, k+1]-1]);
            }
          }
        } else if (link == 3) { // cloglog link
          for (i in 1:ni_agd_arm) {
            if (which_RE[narm_ipd + i])
              eta_agd_arm_RE = eta_agd_arm_noRE[i] + f_delta[which_RE[narm_ipd + i]];
            else
              eta_agd_arm_RE = eta_agd_arm_noRE[i];

            for (k in 1:(agd_arm_ncat[i] - 1)) {
              q_agd_arm_bar[i, k] = inv_cloglog(eta_agd_arm_RE - cc[agd_arm_cat[i, k+1]-1]);
            }
          }
        }

      } else {

        if (class_effects) {
          for (i in 1:ni_agd_arm) {
            if (agd_arm_trt[i] > 1 && which_CE[agd_arm_trt[i] - 1]) {
              eta_agd_arm_noRE[i] += f_class[which_fclass[agd_arm_trt[i] - 1]];
            }
          }
        }

        if (link == 1) { // logit link
          for (i in 1:ni_agd_arm) {
            for (k in 1:(agd_arm_ncat[i] - 1)) {
              q_agd_arm_bar[i, k] = inv_logit(eta_agd_arm_noRE[i] - cc[agd_arm_cat[i, k+1]-1]);
            }
          }
        } else if (link == 2) { // probit link
          for (i in 1:ni_agd_arm) {
            for (k in 1:(agd_arm_ncat[i] - 1)) {
              q_agd_arm_bar[i, k] = Phi(eta_agd_arm_noRE[i] - cc[agd_arm_cat[i, k+1]-1]);
            }
          }
        } else if (link == 3) { // cloglog link
          for (i in 1:ni_agd_arm) {
            for (k in 1:(agd_arm_ncat[i] - 1)) {
              q_agd_arm_bar[i, k] = inv_cloglog(eta_agd_arm_noRE[i] - cc[agd_arm_cat[i, k+1]-1]);
            }
          }
        }
      }
    }

    for (i in 1:ni_agd_arm) {
      // Category 1
      theta_agd_arm_bar[i, 1] = 1 - q_agd_arm_bar[i, 1];

      // Categories 2:(agd_arm_ncat - 1)
      for (k in 2:(agd_arm_ncat[i] - 1))
        theta_agd_arm_bar[i, agd_arm_cat[i, k]] = q_agd_arm_bar[i, k - 1] - q_agd_arm_bar[i, k];

      // Category agd_arm_ncat
      theta_agd_arm_bar[i, agd_arm_cat[i, agd_arm_ncat[i]]] = q_agd_arm_bar[i, agd_arm_ncat[i] - 1];
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

            mu_ful[(c_x+1):(c_x+agd_regression_nx[i])] =  1 - inv_logit(X_agd_regression_int[ (c_x+1):(c_x+agd_regression_nx[i]),] * allbeta) ;

            mu_err[(c_x+1):(c_x+agd_regression_nx[i])] =
              mu_red[ (c_x+1):(c_x+agd_regression_nx[i])] -
              mu_ful[ (c_x+1):(c_x+agd_regression_nx[i])];

            lp_err[(c_x+1):(c_x+agd_regression_nx[i])] =
              eta_red[(c_x+1):(c_x+agd_regression_nx[i])] -
              logit(mu_red[(c_x+1):(c_x+agd_regression_nx[i])] - mu_err[(c_x+1):(c_x+agd_regression_nx[i])]);

            eta_agd_regression[(c_c+1):(c_c+agd_regression_ncoef[i])] =
            (
              allbeta[XI_col_vec[(c_i+1):(c_i+agd_regression_ncoef_inc[i])]] +
              block(agd_regression_OVB_mat_omt[i], 1, 1,agd_regression_ncoef_inc[i] ,agd_regression_ncoef_omt[i] ) *
              allbeta[XO_col_vec[(c_o+1):(c_o+agd_regression_ncoef_omt[i])]] +
              block(agd_regression_OVB_mat_hat[i]  , 1, 1,agd_regression_ncoef_inc[i] ,agd_regression_nx[i] ) *
              lp_err[(c_x+1):(c_x+agd_regression_nx[i])]
            )[agd_regression_est_map_cc] + X_agd_regression_cc[(c_c+1):(c_c+agd_regression_ncoef[i]), ] * cc;

          }else if (link == 2){ // probit link

            mu_ful[(c_x+1):(c_x+agd_regression_nx[i])] =  1 - Phi(X_agd_regression_int[ (c_x+1):(c_x+agd_regression_nx[i]),] * allbeta) ;

            mu_err[(c_x+1):(c_x+agd_regression_nx[i])] =
              mu_red[ (c_x+1):(c_x+agd_regression_nx[i])] -
              mu_ful[ (c_x+1):(c_x+agd_regression_nx[i])];

            lp_err[(c_x+1):(c_x+agd_regression_nx[i])] =
              eta_red[(c_x+1):(c_x+agd_regression_nx[i])] -
              inv_Phi(mu_red[(c_x+1):(c_x+agd_regression_nx[i])] - mu_err[(c_x+1):(c_x+agd_regression_nx[i])]);

            eta_agd_regression[(c_c+1):(c_c+agd_regression_ncoef[i])] =
            (
              allbeta[XI_col_vec[(c_i+1):(c_i+agd_regression_ncoef_inc[i])]] +
              block(agd_regression_OVB_mat_omt[i], 1, 1,agd_regression_ncoef_inc[i] ,agd_regression_ncoef_omt[i] ) *
              allbeta[XO_col_vec[(c_o+1):(c_o+agd_regression_ncoef_omt[i])]] +
              block(agd_regression_OVB_mat_hat[i]  , 1, 1,agd_regression_ncoef_inc[i] ,agd_regression_nx[i] ) *
              lp_err[(c_x+1):(c_x+agd_regression_nx[i])]
            )[agd_regression_est_map_cc] + X_agd_regression_cc[(c_c+1):(c_c+agd_regression_ncoef[i]), ] * cc;

          }else if (link == 3){ // cloglog link

            mu_ful[(c_x+1):(c_x+agd_regression_nx[i])] =  1 - inv_cloglog(X_agd_regression_int[ (c_x+1):(c_x+agd_regression_nx[i]),] * allbeta) ;

            mu_err[(c_x+1):(c_x+agd_regression_nx[i])] =
              mu_red[ (c_x+1):(c_x+agd_regression_nx[i])] -
              mu_ful[ (c_x+1):(c_x+agd_regression_nx[i])];

            lp_err[(c_x+1):(c_x+agd_regression_nx[i])] =
              eta_red[(c_x+1):(c_x+agd_regression_nx[i])] -
              log(-log1m(mu_red[(c_x+1):(c_x+agd_regression_nx[i])] - mu_err[(c_x+1):(c_x+agd_regression_nx[i])]));

            eta_agd_regression[(c_c+1):(c_c+agd_regression_ncoef[i])] =
            (
              allbeta[XI_col_vec[(c_i+1):(c_i+agd_regression_ncoef_inc[i])]] +
              block(agd_regression_OVB_mat_omt[i], 1, 1,agd_regression_ncoef_inc[i] ,agd_regression_ncoef_omt[i] ) *
              allbeta[XO_col_vec[(c_o+1):(c_o+agd_regression_ncoef_omt[i])]] +
              block(agd_regression_OVB_mat_hat[i]  , 1, 1,agd_regression_ncoef_inc[i] ,agd_regression_nx[i] ) *
              lp_err[(c_x+1):(c_x+agd_regression_nx[i])]
            )[agd_regression_est_map_cc] + X_agd_regression_cc[(c_c+1):(c_c+agd_regression_ncoef[i]), ] * cc;

          }
        } else {
            eta_agd_regression[ (c_c+1):(c_c+agd_regression_ncoef[i]) ] = X_agd_regression_no_QR[(c_c+1):(c_c+agd_regression_ncoef[i]), ] * allbeta +
                                                                          X_agd_regression_cc[(c_c+1):(c_c+agd_regression_ncoef[i]), ] * cc;
        }
        c_c += agd_regression_ncoef[i];
        c_i += agd_regression_ncoef_inc[i];
        c_o += agd_regression_ncoef_omt[i];
        c_x += agd_regression_nx[i];
      }
    }else{
      eta_agd_regression = X_agd_regression * beta_tilde + X_agd_regression_cc * cc ;
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
  for (i in 1:ni_ipd)
    ipd_r[i] ~ categorical(theta_ipd[i]);
  // Use ordered_(logistic|probit) for those links?

  // -- AgD likelihood (arm-based) --
  for (i in 1:ni_agd_arm)
    agd_arm_r[i] ~ multinomial(theta_agd_arm_bar[i]);

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

  // -- Priors on cutpoints --
  // Implied improper uniform prior on cutpoints if prior_aux_dist = 0
  //   cc ~ uniform(-inf, inf)
  // Otherwise put priors on the differences between cutpoints. Stan will
  // automatically impose the ordering constraints
  if (prior_aux_dist > 0) {
    vector[ncat - 2] diff_cc;
    for (k in 1:(ncat - 2))
      diff_cc[k] = cc[k + 1] - cc[k];
    prior_select_lp(diff_cc, prior_aux_dist, prior_aux_location, prior_aux_scale, prior_aux_df);
  }
}
generated quantities {
  // Note: fitted values and theta_bar_cum_agd_arm will be 0 for missing categories
  array[ni_ipd] vector[ncat] fitted_ipd;
  array[ni_agd_arm] vector[ncat] fitted_agd_arm;
  matrix[ni_agd_arm * n_int_thin, ncat] theta_bar_cum_agd_arm = rep_matrix(0, ni_agd_arm * n_int_thin, ncat);

#include /include/generated_quantities_common.stan

  // IPD log likelihood and residual deviance
  for (i in 1:ni_ipd) {
    log_lik[i] = categorical_lpmf(ipd_r[i] | theta_ipd[i]);
    resdev[i] = -2 * log_lik[i];
    fitted_ipd[i] = theta_ipd[i];
  }

  // AgD (arm-based) log likelihood and residual deviance
  {
    array[ni_agd_arm] vector[ncat] dv;
    for (i in 1:ni_agd_arm) {
      log_lik[ni_ipd + i] = multinomial_lpmf(agd_arm_r[i] | theta_agd_arm_bar[i]);
      fitted_agd_arm[i] = agd_arm_n[i] * theta_agd_arm_bar[i];

      for (k in 1:agd_arm_ncat[i]) {
        // Multinomial residual deviance
        dv[i, k] = agd_arm_r[i, agd_arm_cat[i, k]] == 0 ? 0 : lmultiply(agd_arm_r[i, agd_arm_cat[i, k]], agd_arm_r[i, agd_arm_cat[i, k]] / fitted_agd_arm[i, agd_arm_cat[i, k]]);
      }
      resdev[ni_ipd + i] = 2 * sum(dv[i, 1:agd_arm_ncat[i]]);
    }
  }

  // Cumulative integration
  for (i in 1:ni_agd_arm) {
    for (j in 1:n_int_thin) {
      vector[ncat - 1] q_agd_arm_bar_cum;
      for (k in 1:(agd_arm_ncat[i] - 1)) q_agd_arm_bar_cum[k] = mean(theta_agd_arm_ii[(1 + (i - 1)*nint):((i - 1)*nint + j*int_thin), k]);

      // Category 1
      theta_bar_cum_agd_arm[(i - 1)*n_int_thin + j, 1] = 1 - q_agd_arm_bar_cum[1];

      // Categories 2:(agd_arm_ncat - 1)
      for (k in 2:(agd_arm_ncat[i] - 1))
        theta_bar_cum_agd_arm[(i - 1)*n_int_thin + j, agd_arm_cat[i, k]] = q_agd_arm_bar_cum[k - 1] - q_agd_arm_bar_cum[k];

      // Category agd_arm_ncat
      theta_bar_cum_agd_arm[(i - 1)*n_int_thin + j, agd_arm_cat[i, agd_arm_ncat[i]]] = q_agd_arm_bar_cum[agd_arm_ncat[i] - 1];
    }
  }
}
