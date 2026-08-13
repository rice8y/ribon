/* Test-only ViennaRNA API adapter; contains no prediction implementation. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <ViennaRNA/eval.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/structures/centroid.h>
#include <ViennaRNA/structures/mea.h>
#include <ViennaRNA/structures/problist.h>
#include <ViennaRNA/probabilities/structures.h>
#include <ViennaRNA/utils/basic.h>

int main(int argc, char **argv) {
  if (argc < 2 || argc > 4) {
    fputs("usage: vienna_fold_harness RNA_SEQUENCE [DANGLES] [STRUCTURE|--verbose]\n", stderr);
    return 2;
  }
  const char *sequence = argv[1];
  int dangles = (argc >= 3) ? atoi(argv[2]) : 2;
  if (dangles < 0 || dangles > 3) {
    fputs("DANGLES must be 0, 1, 2, or 3\n", stderr);
    return 2;
  }
  size_t length = strlen(sequence);
  char *mfe_structure = vrna_alloc(length + 1);
  char *probability_string = vrna_alloc(length + 1);
  vrna_md_t md;
  vrna_md_set_default(&md);
  md.dangles = dangles;
  md.uniq_ML = 1;
  vrna_fold_compound_t *fc = vrna_fold_compound(sequence, &md, VRNA_OPTION_DEFAULT);
  double mfe = vrna_mfe(fc, mfe_structure);
  const char *evaluated_structure =
      (argc == 4 && strncmp(argv[3], "--", 2) != 0) ? argv[3] : mfe_structure;
  double evaluated_mfe = vrna_eval_structure(fc, mfe_structure);
  double evaluated_supplied = vrna_eval_structure(fc, evaluated_structure);
  short *evaluated_pt = vrna_ptable(evaluated_structure);
  int evaluated_exterior = vrna_eval_loop_pt(fc, 0, evaluated_pt);
  if (argc == 4 && strcmp(argv[3], "--verbose") == 0)
    vrna_eval_structure_verbose(fc, mfe_structure, stderr);
  vrna_exp_params_rescale(fc, &mfe);
  double ensemble_energy = vrna_pf(fc, probability_string);
  double centroid_distance = 0.0;
  char *centroid = vrna_centroid(fc, &centroid_distance);
  float mea_score = 0.0f;
  char *mea = vrna_MEA(fc, 1.0, &mea_score);
  vrna_ep_t *probabilities = vrna_plist_from_probs(fc, 1e-12);
  double mean_bp_distance = vrna_mean_bp_distance(fc);
  double *positional_entropy = vrna_positional_entropy(fc);
  double ensemble_defect = vrna_ensemble_defect(fc, mfe_structure);

  printf("{\"sequence\":\"%s\",\"dangles\":%d,\"uniq_ML\":true,", sequence, dangles);
  printf("\"mfe_structure\":\"%s\",\"mfe_energy\":%.9g,", mfe_structure, mfe);
  printf("\"evaluated_mfe_energy\":%.9g,", evaluated_mfe);
  printf("\"evaluated_structure\":\"%s\",\"evaluated_structure_energy\":%.9g,",
         evaluated_structure, evaluated_supplied);
  printf("\"evaluated_exterior_energy\":%.9g,", evaluated_exterior / 100.0);
  printf("\"evaluated_loop_energies\":[");
  int first_loop = 1;
  for (size_t position = 1; position <= length; ++position) {
    if (evaluated_pt[position] <= (short)position) continue;
    int loop_energy = vrna_eval_loop_pt(fc, (int)position, evaluated_pt);
    if (!first_loop) putchar(',');
    first_loop = 0;
    printf("{\"i\":%zu,\"j\":%d,\"energy\":%.9g}",
           position, evaluated_pt[position], loop_energy / 100.0);
  }
  printf("],");
  printf("\"ensemble_free_energy\":%.9g,", ensemble_energy);
  printf("\"centroid_structure\":\"%s\",\"centroid_distance\":%.9g,", centroid, centroid_distance);
  printf("\"mea_structure\":\"%s\",\"mea_score\":%.9g,", mea, mea_score);
  printf("\"mean_base_pair_distance\":%.9g,", mean_bp_distance);
  printf("\"mfe_ensemble_defect\":%.9g,", ensemble_defect);
  printf("\"positional_entropy_bits\":[");
  for (size_t position = 1; position <= length; ++position) {
    if (position > 1) putchar(',');
    printf("%.9g", positional_entropy[position]);
  }
  printf("],");
  printf("\"pair_probabilities\":[");
  int first = 1;
  for (vrna_ep_t *entry = probabilities; entry && entry->i != 0; ++entry) {
    if (entry->type != VRNA_PLIST_TYPE_BASEPAIR) continue;
    if (!first) putchar(',');
    first = 0;
    printf("{\"i\":%d,\"j\":%d,\"p\":%.9g}", entry->i, entry->j, entry->p);
  }
  puts("]}");

  free(probabilities);
  free(positional_entropy);
  free(evaluated_pt);
  free(mea);
  free(centroid);
  free(probability_string);
  free(mfe_structure);
  vrna_fold_compound_free(fc);
  return 0;
}
