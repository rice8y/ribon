/* Test-only adapter for the ViennaRNA RNAduplex MFE implementation. */
#include <stdio.h>
#include <stdlib.h>

#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/duplex.h>
#include <ViennaRNA/model.h>

int main(int argc, char **argv) {
  if (argc < 3 || argc > 4) {
    fputs("usage: vienna_duplex_harness SEQUENCE_A SEQUENCE_B [SALT_MOLAR]\n", stderr);
    return 2;
  }

  vrna_md_defaults_dangles(0);
  if (argc == 4)
    vrna_md_defaults_salt(strtod(argv[3], NULL));

  duplexT result = duplexfold(argv[1], argv[2]);
  printf("{\"sequence_a\":\"%s\",\"sequence_b\":\"%s\",", argv[1], argv[2]);
  printf("\"salt_molar\":%.12g,\"structure\":\"%s\",", argc == 4 ? strtod(argv[3], NULL) : 1.021, result.structure);
  printf("\"energy\":%.12g,\"i\":%d,\"j\":%d,\"suboptimal\":[", result.energy, result.i, result.j);
  duplexT *suboptimal = duplex_subopt(argv[1], argv[2], 1000, 0);
  for (int index = 0; suboptimal[index].i != 0; ++index) {
    if (index) putchar(',');
    printf("{\"structure\":\"%s\",\"energy\":%.12g}", suboptimal[index].structure, suboptimal[index].energy);
    free(suboptimal[index].structure);
  }
  puts("]}");
  free(suboptimal);
  free(result.structure);
  return 0;
}
