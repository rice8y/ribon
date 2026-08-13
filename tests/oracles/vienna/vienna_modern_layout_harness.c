/* Test-only adapter for ViennaRNA RNAturtle/RNApuzzler coordinate APIs. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ViennaRNA/plotting/RNApuzzler/RNApuzzler.h>
#include <ViennaRNA/plotting/RNApuzzler/RNAturtle.h>

int main(int argc, char **argv) {
  if (argc != 3 || (strcmp(argv[1], "turtle") != 0 && strcmp(argv[1], "puzzler") != 0)) {
    fputs("usage: vienna_modern_layout_harness turtle|puzzler DOT_BRACKET\n", stderr);
    return 2;
  }
  float *x = NULL;
  float *y = NULL;
  double *arcs = NULL;
  int length = strcmp(argv[1], "turtle") == 0
                   ? vrna_plot_coords_turtle(argv[2], &x, &y, &arcs)
                   : vrna_plot_coords_puzzler(argv[2], &x, &y, &arcs, NULL);
  if (length <= 0) return 3;
  int finite = 1;
  for (int i = 0; i < length; ++i)
    if (!isfinite(x[i]) || !isfinite(y[i])) finite = 0;
  printf("{\"method\":\"%s\",\"structure\":\"%s\",\"finite\":%s,\"points\":[",
         argv[1], argv[2], finite ? "true" : "false");
  for (int i = 0; i < length; ++i) {
    if (i) putchar(',');
    if (isfinite(x[i]) && isfinite(y[i]))
      printf("{\"x\":%.9g,\"y\":%.9g}", x[i], y[i]);
    else
      fputs("{\"x\":null,\"y\":null}", stdout);
  }
  puts("]}");
  free(arcs);
  free(x);
  free(y);
  return 0;
}
