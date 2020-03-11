// Function for selecting given prior distribution, provides sampling statement ~prior_select(...)
void prior_select_lp(vector y, int dist, real location, real scale, real df) {
  if (dist == 1) { // Normal
    y ~ normal(location, scale);
  } else if (dist == 2) { // Cauchy
    y ~ cauchy(location, scale);
  } else if (dist == 3) { // Student t
    y ~ student_t(df, location, scale);
  } else if (dist == 4) { // Exponential
    y ~ exponential(1/scale);
  } else {
    reject("Not a supported prior dist.");
  }
  return;
}
