// Function for selecting given prior distribution, provides sampling statement ~prior_select(...)
void prior_select_lp(vector y, int type, real location, real scale, real df) {
  if (type == 1) { // Normal
    y ~ normal(location, scale);
  } else if (type == 2) { // Cauchy
    y ~ cauchy(location, scale);
  } else if (type == 3) { // Student t
    y ~ student_t(df, location, scale);
  } else if (type == 4) { // Exponential
    y ~ exponential(1/scale);
  } else {
    reject("Not a supported prior type.");
  }
  return;
}
