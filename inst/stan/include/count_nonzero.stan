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
