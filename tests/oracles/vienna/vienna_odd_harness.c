/* Test-only ViennaRNA API adapter; contains no prediction implementation. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <ViennaRNA/eval.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/utils/basic.h>

int main(int argc, char **argv) {
  if (argc < 3 || argc > 4) {
    fputs("usage: vienna_odd_harness SEQUENCE DANGLES [STRUCTURE]\n", stderr);
    return 2;
  }
  const char *sequence = argv[1];
  int dangles = atoi(argv[2]);
  if (dangles != 1 && dangles != 3) {
    fputs("DANGLES must be 1 or 3\n", stderr);
    return 2;
  }

  size_t length = strlen(sequence);
  char *mfe_structure = vrna_alloc(length + 1);
  vrna_md_t md;
  vrna_md_set_default(&md);
  md.dangles = dangles;
  md.uniq_ML = 1;
  vrna_fold_compound_t *fc =
      vrna_fold_compound(sequence, &md, VRNA_OPTION_DEFAULT);
  double mfe = vrna_mfe(fc, mfe_structure);
  const char *evaluated_structure = argc == 4 ? argv[3] : mfe_structure;
  double evaluated = vrna_eval_structure(fc, evaluated_structure);
  short *pt = vrna_ptable(evaluated_structure);

  printf("{\"sequence\":\"%s\",\"dangles\":%d,", sequence, dangles);
  printf("\"mfe_structure\":\"%s\",\"mfe_energy\":%.9g,", mfe_structure,
         mfe);
  printf("\"evaluated_structure\":\"%s\",\"evaluated_energy\":%.9g,",
         evaluated_structure, evaluated);
  printf("\"loop_energies\":[{\"i\":0,\"j\":0,\"energy\":%.9g}",
         vrna_eval_loop_pt(fc, 0, pt) / 100.0);
  for (size_t i = 1; i <= length; ++i) {
    if (pt[i] <= (short)i)
      continue;
    printf(",{\"i\":%zu,\"j\":%d,\"energy\":%.9g}", i, pt[i],
           vrna_eval_loop_pt(fc, (int)i, pt) / 100.0);
  }
  puts("]}");

  free(pt);
  free(mfe_structure);
  vrna_fold_compound_free(fc);
  return 0;
}
