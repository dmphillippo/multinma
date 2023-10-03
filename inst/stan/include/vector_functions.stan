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
