// -- Vector functions --

// Which entries are equal to a given value
int[] which(int[] x, int y) {
  int n = num_elements(x);
  int w[n]; // Over-allocate w and then truncate later
  int c = 1;
  for (i in 1:n) {
    if (x[i] == y) {
      w[c] = i;
      c += 1;
    }
  }
  return w[1:(c-1)];
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

// pow
vector pow_vec (vector x, vector y) {
  // (Not needed past Stan 2.26)
  int n = num_elements(x);
  vector[n] out;
  for (i in 1:n) out[i] = pow(x[i], y[i]);
  return out;
}

vector pow_vec2 (vector x, real y) {
  // (Not needed past Stan 2.26)
  int n = num_elements(x);
  vector[n] out;
  for (i in 1:n) out[i] = pow(x[i], y);
  return out;
}

// lmultiply
vector lmultiply_vec (vector x, vector y) {
  int n = num_elements(x);
  vector[n] out;
  for (i in 1:n) out[i] = lmultiply(x[i], y[i]);
  return out;
}
