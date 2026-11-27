/* Recupera los polinomios del gadget (tg_0, tg_1) usando la implementacion
 * ORIGINAL de D=64. Se compila contra `auditor/` sin tocar nada y escribe por
 * stdout D lineas "poly1 poly2" por cada entrada de las tablas.
 *
 * gentdrow.c se incluye para acceder a las tablas static const; por eso el
 * build debe excluir gentdrow.c de la libreria (SKIP_EXTRA=gentdrow.c).
 */
#include <stdint.h>
#include <stdio.h>

#include "param.h"
#include "poly_qbig.h"

#include "gentdrow.c"

int main(void) {
  POLY_QBIG e;
  uint64_t k, j;

  for (k = 0; k < (uint64_t)(N_SPENT - 1); k++) {
    e = tg_0[k];
    intt_qbig(&e);
    for (j = 0; j < D; j++)
      printf("%llu %llu\n", (unsigned long long)e.poly1[j],
             (unsigned long long)e.poly2[j]);
  }

  for (k = 0; k < (uint64_t)(S * R); k++) {
    e = tg_1[k];
    intt_qbig(&e);
    for (j = 0; j < D; j++)
      printf("%llu %llu\n", (unsigned long long)e.poly1[j],
             (unsigned long long)e.poly2[j]);
  }

  fprintf(stderr, "recuperadas %llu entradas (tg_0=%llu, tg_1=%llu)\n",
          (unsigned long long)((N_SPENT - 1) + S * R),
          (unsigned long long)(N_SPENT - 1), (unsigned long long)(S * R));
  return 0;
}
