/* Test-only adapter for ViennaRNA's circular MFE and partition APIs. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/structures/problist.h>
#include <ViennaRNA/utils/basic.h>

int main(int argc, char **argv) {
  if (argc != 2) {
    fputs("usage: vienna_circular_harness RNA_SEQUENCE\n", stderr);
    return 2;
  }
  size_t length = strlen(argv[1]);
  char *structure = vrna_alloc(length + 1);
  vrna_md_t md;
  vrna_md_set_default(&md);
  md.circ = 1;
  md.dangles = 0;
  md.uniq_ML = 1;
  vrna_fold_compound_t *fc = vrna_fold_compound(argv[1], &md, VRNA_OPTION_DEFAULT);
  double mfe = vrna_mfe(fc, structure);
  vrna_exp_params_rescale(fc, &mfe);
  double ensemble = vrna_pf(fc, NULL);
  vrna_ep_t *probabilities = vrna_plist_from_probs(fc, 1e-12);
  printf("{\"sequence\":\"%s\",\"structure\":\"%s\",", argv[1], structure);
  printf("\"mfe_energy\":%.12g,\"ensemble_free_energy\":%.12g,", mfe, ensemble);
  printf("\"pair_probabilities\":[");
  int first = 1;
  for (vrna_ep_t *entry = probabilities; entry && entry->i != 0; ++entry) {
    if (entry->type != VRNA_PLIST_TYPE_BASEPAIR) continue;
    if (!first) putchar(',');
    first = 0;
    printf("{\"i\":%d,\"j\":%d,\"p\":%.12g}", entry->i, entry->j, entry->p);
  }
  puts("]}");
  free(probabilities);
  free(structure);
  vrna_fold_compound_free(fc);
  return 0;
}
