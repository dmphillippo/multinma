// Function for selecting given prior distribution, provides sampling statement ~prior_select(...)
void prior_select_lp(vector y, int dist, real location, real scale, real df) {
  if (dist == 0) { // Implicit flat prior

  } else if (dist == 1) { // Normal
    y ~ normal(location, scale);
  } else if (dist == 2) { // Cauchy
    y ~ cauchy(location, scale);
  } else if (dist == 3) { // Student t
    y ~ student_t(df, location, scale);
  } else if (dist == 4) { // Exponential
    y ~ exponential(1/scale);
  } else if (dist == 5) { // log-Normal
    y ~ lognormal(location, scale);
  } else if (dist == 6) { // log-Student t
    log(y) ~ student_t(df, location, scale);
    target += -log(y);
  } else {
    reject("Not a supported prior dist.");
  }
  return;
}

void prior_select2_lp(real y, int dist, real location, real scale, real df) {
  if (dist == 0) { // Implicit flat prior

  } else if (dist == 1) { // Normal
    y ~ normal(location, scale);
  } else if (dist == 2) { // Cauchy
    y ~ cauchy(location, scale);
  } else if (dist == 3) { // Student t
    y ~ student_t(df, location, scale);
  } else if (dist == 4) { // Exponential
    y ~ exponential(1/scale);
  } else if (dist == 5) { // log-Normal
    y ~ lognormal(location, scale);
  } else if (dist == 6) { // log-Student t
    log(y) ~ student_t(df, location, scale);
    target += -log(y);
  } else {
    reject("Not a supported prior dist.");
  }
  return;
}

