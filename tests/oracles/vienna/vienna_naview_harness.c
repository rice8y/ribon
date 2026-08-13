/*
 * Test-only adapter for an independently built ViennaRNA naview object.
 *
 * This file contains no ViennaRNA implementation code. It supplies the small
 * allocation/logging surface required to execute vrna_plot_coords_naview_pt()
 * and emits reference coordinates for comparison tests.
 */
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <ViennaRNA/utils/log.h>
#include <ViennaRNA/plotting/naview/naview.h>

void *vrna_alloc(size_t size) {
  void *memory = calloc(1, size);
  if (!memory) {
    fputs("allocation failed\n", stderr);
    exit(2);
  }
  return memory;
}

void vrna_log(vrna_log_levels_e level,
              const char *file_name,
              int line_number,
              const char *format_string,
              ...) {
  (void)level;
  (void)file_name;
  (void)line_number;
  (void)format_string;
}

short *vrna_ptable(const char *structure) {
  (void)structure;
  return NULL;
}

static short *pair_table(const char *structure) {
  size_t length = strlen(structure);
  short *table = calloc(length + 1, sizeof(short));
  int *stack = calloc(length, sizeof(int));
  int top = 0;
  if (!table || !stack || length > 32767) return NULL;
  table[0] = (short)length;
  for (size_t i = 0; i < length; ++i) {
    if (structure[i] == '(') {
      stack[top++] = (int)i + 1;
    } else if (structure[i] == ')') {
      if (top == 0) return NULL;
      int left = stack[--top];
      table[left] = (short)(i + 1);
      table[i + 1] = (short)left;
    } else if (structure[i] != '.') {
      return NULL;
    }
  }
  free(stack);
  if (top != 0) {
    free(table);
    return NULL;
  }
  return table;
}

int main(int argc, char **argv) {
  if (argc != 2) {
    fputs("usage: vienna_naview_harness DOT_BRACKET\n", stderr);
    return 2;
  }
  short *table = pair_table(argv[1]);
  float *x = NULL;
  float *y = NULL;
  if (!table) {
    fputs("invalid dot-bracket string\n", stderr);
    return 2;
  }
  int length = vrna_plot_coords_naview_pt(table, &x, &y);
  if (length <= 0) return 3;
  printf("{\"structure\":\"%s\",\"points\":[", argv[1]);
  for (int i = 0; i < length; ++i) {
    if (i) putchar(',');
    printf("{\"x\":%.9g,\"y\":%.9g}", x[i], y[i]);
  }
  puts("]}");
  free(table);
  free(x);
  free(y);
  return 0;
}

