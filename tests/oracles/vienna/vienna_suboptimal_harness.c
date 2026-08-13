/* Test-only adapter for ViennaRNA's Wuchty suboptimal enumeration. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/subopt/basic.h>
#include <ViennaRNA/subopt/wuchty.h>

int main(int argc, char **argv) {
  if (argc < 3 || argc > 4) {
    fputs("usage: vienna_suboptimal_harness SEQUENCE BAND_KCAL [DANGLES]\n", stderr);
    return 2;
  }
  const char *sequence = argv[1];
  const double band = strtod(argv[2], NULL);
  const int dangles = argc == 4 ? atoi(argv[3]) : 2;
  vrna_md_t md;
  vrna_md_set_default(&md);
  md.dangles = dangles;
  md.uniq_ML = 1;
  vrna_fold_compound_t *fc = vrna_fold_compound(sequence, &md, VRNA_OPTION_MFE);
  vrna_subopt_solution_t *solutions = vrna_subopt(fc, (int)(band * 100.0 + 0.5), 1, NULL);

  printf("{\"sequence\":\"%s\",\"energy_band\":%.12g,\"dangles\":%d,\"structures\":[", sequence, band, dangles);
  for (int index = 0; solutions[index].structure; ++index) {
    if (index) putchar(',');
    printf("{\"structure\":\"%s\",\"energy\":%.12g}", solutions[index].structure, solutions[index].energy);
    free(solutions[index].structure);
  }
  puts("]}");
  free(solutions);
  vrna_fold_compound_free(fc);
  return 0;
}
