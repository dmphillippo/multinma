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
