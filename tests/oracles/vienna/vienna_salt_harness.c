/* Test-only adapter for ViennaRNA non-default monovalent-salt parameters. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <ViennaRNA/eval.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/structures/problist.h>
#include <ViennaRNA/utils/basic.h>

int main(int argc, char **argv) {
  if (argc < 3 || argc > 5) {
    fputs("usage: vienna_salt_harness SEQUENCE SALT_MOLAR [DANGLES] [STRUCTURE]\n", stderr);
    return 2;
  }
  const char *sequence = argv[1];
  const double salt = strtod(argv[2], NULL);
  const int dangles = argc >= 4 ? atoi(argv[3]) : 2;
  const size_t length = strlen(sequence);
  char *mfe_structure = vrna_alloc(length + 1);
  char *pf_structure = vrna_alloc(length + 1);
  vrna_md_t md;
  vrna_md_set_default(&md);
  md.salt = salt;
  md.dangles = dangles;
  md.uniq_ML = 1;
  vrna_fold_compound_t *fc = vrna_fold_compound(sequence, &md, VRNA_OPTION_DEFAULT);
  double mfe = vrna_mfe(fc, mfe_structure);
  const char *evaluated = argc == 5 ? argv[4] : mfe_structure;
  const double evaluated_energy = vrna_eval_structure(fc, evaluated);
  vrna_exp_params_rescale(fc, &mfe);
  const double ensemble = vrna_pf(fc, pf_structure);
  vrna_ep_t *probabilities = vrna_plist_from_probs(fc, 1e-10);

  printf("{\"sequence\":\"%s\",\"salt_molar\":%.12g,\"dangles\":%d,", sequence, salt, dangles);
  printf("\"mfe_structure\":\"%s\",\"mfe_energy\":%.12g,", mfe_structure, mfe);
  printf("\"evaluated_structure\":\"%s\",\"evaluated_energy\":%.12g,", evaluated, evaluated_energy);
  printf("\"ensemble_free_energy\":%.12g,\"pair_probabilities\":[", ensemble);
  int first = 1;
  for (vrna_ep_t *entry = probabilities; entry && entry->i != 0; ++entry) {
    if (entry->type != VRNA_PLIST_TYPE_BASEPAIR) continue;
    if (!first) putchar(',');
    first = 0;
    printf("{\"i\":%d,\"j\":%d,\"p\":%.12g}", entry->i, entry->j, entry->p);
  }
  puts("]}");

  free(probabilities);
  free(pf_structure);
  free(mfe_structure);
  vrna_fold_compound_free(fc);
  return 0;
}
