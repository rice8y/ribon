/* Test-only adapter for ViennaRNA's public G-quadruplex energy API. */
#include <stdio.h>
#include <stdlib.h>

#include <ViennaRNA/eval/gquad.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/params/basic.h>

int main(int argc, char **argv) {
  if (argc != 6) {
    fputs("usage: vienna_gquad_harness STACK LINKER1 LINKER2 LINKER3 TEMPERATURE_C\n", stderr);
    return 2;
  }
  unsigned int stack = (unsigned int)strtoul(argv[1], NULL, 10);
  unsigned int linkers[3] = {
      (unsigned int)strtoul(argv[2], NULL, 10),
      (unsigned int)strtoul(argv[3], NULL, 10),
      (unsigned int)strtoul(argv[4], NULL, 10),
  };
  vrna_md_t md;
  vrna_md_set_default(&md);
  md.temperature = strtod(argv[5], NULL);
  md.gquad = 1;
  vrna_param_t *parameters = vrna_params(&md);
  int energy = vrna_E_gquad(stack, linkers, parameters);
  printf("{\"stack_size\":%u,\"linkers\":[%u,%u,%u],\"temperature_celsius\":%.12g,\"energy_kcal_mol\":%.12g}\n",
         stack, linkers[0], linkers[1], linkers[2], md.temperature, energy / 100.0);
  free(parameters);
  return 0;
}
