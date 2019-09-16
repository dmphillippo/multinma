// Function for selecting given prior distribution, provides sampling statement ~prior_select(...)
real prior_select_lpdf(real y, int type, real location, real scale, real df) {
  if (type == 1) { // Normal
    return normal_lpdf(y | location, scale);
  } else if (type == 2) { // Cauchy
    return cauchy_lpdf(y | location, scale);
  } else if (type == 3) { // Student t
    return student_t_lpdf(y | df, location, scale);
  } else if (type == 4) { // Exponential
    return exponential_lpdf(y | 1/scale);
  }
}
