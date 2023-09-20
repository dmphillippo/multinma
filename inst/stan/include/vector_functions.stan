// -- Vector functions --

// Which entries are equal to a given value
array[] int which(array[] int x, int y) {
  int n = num_elements(x);
  array[n] int w; // Over-allocate w and then truncate later
  int c = 1;
  for (i in 1:n) {
    if (x[i] == y) {
      w[c] = i;
      c += 1;
    }
  }
  return w[1:(c-1)];
}

// Which entry is first equal to a given value
int first(array[] int x, int y) {
  int w = 0;
  int i = 1;
  while (i <= num_elements(x) && w == 0) {
    if (x[i] == y) w = i;
    i += 1;
  }
  return w;
}

// How many entries are equal to a given value
int nwhich(array[] int x, int y) {
  int w = 0;
  for (i in 1:num_elements(x)) {
    if (x[i] == y) w += 1;
  }
  return w;
}

array[] int nwhich_all(array[] int x, int max_id) {
  array[max_id] int w = rep_array(0, max_id);
  for (i in 1:num_elements(x)) {
    if (x[i]) w[x[i]] += 1;
  }
  return w;
}

// Which entries are greater than 0
array[] int which_gt0(vector x) {
  int n = num_elements(x);
  array[n] int w; // Over-allocate w and then truncate later
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
