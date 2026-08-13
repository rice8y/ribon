/* Test-only adapter for ViennaRNA's RNAplfold-compatible public wrappers. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <ViennaRNA/model.h>
#include <ViennaRNA/partfunc/local.h>

int main(int argc, char **argv) {
  if (argc != 5) {
    fputs("usage: vienna_local_harness SEQUENCE WINDOW SPAN UNPAIRED\n", stderr);
    return 2;
  }
  const char *sequence = argv[1];
  int window = atoi(argv[2]);
  int span = atoi(argv[3]);
  int unpaired = atoi(argv[4]);
  int length = (int)strlen(sequence);
  vrna_md_defaults_dangles(0);
  vrna_ep_t *pairs = vrna_pfl_fold(sequence, window, span, 1e-12f);
  double **up = unpaired > 0 ? vrna_pfl_fold_up(sequence, unpaired, window, span) : NULL;
  printf("{\"sequence\":\"%s\",\"window\":%d,\"span\":%d,", sequence, window, span);
  printf("\"pair_probabilities\":[");
  int first = 1;
  for (vrna_ep_t *entry = pairs; entry && entry->i != 0; ++entry) {
    if (entry->type != VRNA_PLIST_TYPE_BASEPAIR) continue;
    if (!first) putchar(',');
    first = 0;
    printf("{\"i\":%d,\"j\":%d,\"p\":%.12g}", entry->i, entry->j, entry->p);
  }
  printf("],\"accessibility\":[");
  first = 1;
  if (up) {
    /* The public matrix stores an interval of length u at its
       right endpoint i (RNAplfold's callback API uses the same coordinates). */
    for (int i = 1; i <= length; ++i) {
      for (int u = 1; u <= unpaired && u <= i; ++u) {
        if (!first) putchar(',');
        first = 0;
        printf("{\"from\":%d,\"length\":%d,\"p\":%.12g}", i - u + 1, u, up[i][u]);
      }
    }
  }
  puts("]}");
  if (up) {
    for (int i = 1; i <= length; ++i) free(up[i]);
    free(up);
  }
  free(pairs);
  return 0;
}
