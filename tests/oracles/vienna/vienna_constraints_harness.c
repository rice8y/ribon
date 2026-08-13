/* Test-only ViennaRNA constraint/probing adapter; no prediction code copied. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <ViennaRNA/constraints/hard.h>
#include <ViennaRNA/constraints/soft.h>
#include <ViennaRNA/eval.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/probing/SHAPE.h>
#include <ViennaRNA/structures/problist.h>
#include <ViennaRNA/utils/basic.h>

static double *parse_reactivities(const char *csv, size_t length) {
  double *values = vrna_alloc((length + 1) * sizeof(double));
  char *copy = strdup(csv);
  char *token = strtok(copy, ",");
  size_t position = 1;
  while (token && position <= length) {
    values[position++] = strcmp(token, "null") == 0 ? -999.0 : atof(token);
    token = strtok(NULL, ",");
  }
  free(copy);
  if (position != length + 1 || token) {
    free(values);
    return NULL;
  }
  return values;
}

int main(int argc, char **argv) {
  if (argc < 8) {
    fputs("usage: vienna_constraints_harness SEQ DANGLES HC NOLP NOGU MAXSPAN MODE [MODE_ARGS]\n",
          stderr);
    return 2;
  }
  const char *sequence = argv[1];
  int dangles = atoi(argv[2]);
  const char *hard = argv[3];
  const char *mode = argv[7];
  size_t length = strlen(sequence);

  vrna_md_t md;
  vrna_md_set_default(&md);
  md.dangles = dangles;
  md.uniq_ML = 1;
  md.noLP = atoi(argv[4]);
  md.noGU = atoi(argv[5]);
  md.max_bp_span = atoi(argv[6]);
  vrna_fold_compound_t *fc =
      vrna_fold_compound(sequence, &md, VRNA_OPTION_DEFAULT);
  if (strcmp(hard, "-") != 0 &&
      !vrna_hc_add_from_db(fc, hard,
                           VRNA_CONSTRAINT_DB_DEFAULT |
                               VRNA_CONSTRAINT_DB_ENFORCE_BP)) {
    fputs("failed to install hard constraints\n", stderr);
    return 2;
  }

  if (strcmp(mode, "none") == 0) {
    /* no-op */
  } else if (strcmp(mode, "up") == 0 && argc == 10) {
    vrna_sc_add_up(fc, (unsigned int)atoi(argv[8]), atof(argv[9]),
                   VRNA_OPTION_DEFAULT);
  } else if (strcmp(mode, "pair") == 0 && argc == 11) {
    vrna_sc_add_bp(fc, (unsigned int)atoi(argv[8]),
                   (unsigned int)atoi(argv[9]), atof(argv[10]),
                   VRNA_OPTION_DEFAULT);
  } else if (strcmp(mode, "deigan") == 0 && argc == 11) {
    double *values = parse_reactivities(argv[10], length);
    if (!values || !vrna_sc_add_SHAPE_deigan(fc, values, atof(argv[8]),
                                              atof(argv[9]),
                                              VRNA_OPTION_DEFAULT)) {
      fputs("failed to install Deigan constraints\n", stderr);
      return 2;
    }
    free(values);
  } else if (strcmp(mode, "zarringhalam") == 0 && argc == 12) {
    double *values = parse_reactivities(argv[11], length);
    if (!values || !vrna_sc_add_SHAPE_zarringhalam(
                       fc, values, atof(argv[8]), atof(argv[10]), argv[9],
                       VRNA_OPTION_DEFAULT)) {
      fputs("failed to install Zarringhalam constraints\n", stderr);
      return 2;
    }
    free(values);
  } else {
    fputs("invalid soft-constraint mode or arguments\n", stderr);
    return 2;
  }

  char *structure = vrna_alloc(length + 1);
  double mfe = vrna_mfe(fc, structure);
  double evaluated = vrna_eval_structure(fc, structure);
  vrna_exp_params_rescale(fc, &mfe);
  char *pf_structure = vrna_alloc(length + 1);
  double ensemble = vrna_pf(fc, pf_structure);
  vrna_ep_t *probabilities = vrna_plist_from_probs(fc, 1e-12);

  printf("{\"structure\":\"%s\",\"mfe\":%.9g,\"evaluated\":%.9g,",
         structure, mfe, evaluated);
  printf("\"ensemble\":%.9g,\"pairs\":[", ensemble);
  int first = 1;
  for (vrna_ep_t *entry = probabilities; entry && entry->i != 0; ++entry) {
    if (entry->type != VRNA_PLIST_TYPE_BASEPAIR)
      continue;
    if (!first)
      putchar(',');
    first = 0;
    printf("{\"i\":%d,\"j\":%d,\"p\":%.9g}", entry->i, entry->j,
           entry->p);
  }
  puts("]}");

  free(probabilities);
  free(pf_structure);
  free(structure);
  vrna_fold_compound_free(fc);
  return 0;
}
