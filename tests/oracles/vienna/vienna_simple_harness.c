/* Test-only adapter for an independently built ViennaRNA library. */
#include <stdio.h>
#include <stdlib.h>

#include <ViennaRNA/plotting/layouts.h>

int main(int argc, char **argv) {
  if (argc != 2) {
    fputs("usage: vienna_simple_harness DOT_BRACKET\n", stderr);
    return 2;
  }
  float *x = NULL;
  float *y = NULL;
  int length = vrna_plot_coords_simple(argv[1], &x, &y);
  if (length <= 0 || !x || !y) return 3;
  printf("{\"structure\":\"%s\",\"points\":[", argv[1]);
  for (int i = 0; i < length; ++i) {
    if (i) putchar(',');
    printf("[%.9g,%.9g]", x[i], y[i]);
  }
  puts("]}");
  free(x);
  free(y);
  return 0;
}
