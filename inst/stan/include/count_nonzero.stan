// Count number of non-zero entries in a matrix
int count_nonzero(matrix m) {
  int I = rows(m);
  int J = cols(m);
  int c = 0;
  for (j in 1:I) {
    for (i in 1:J) {
      if (m[i,j] != 0) c += 1;
    }
  }
  return c;
}

// Which entries are greater than 0
array[] int which_gt0a(array[] int x) {
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
