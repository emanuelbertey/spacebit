/* Prueba REAL de Au con dimension de anillo mayor: D=256 (en vez de 64).
 * Ronda completa: keygen + mint + spend + verify. NO toca fuentes originales.
 * Se compila con -DD=256. No bloquea (sin getchar).
 * Devuelve 0 si verifica, 1 si falla.
 * NOTA: el esquema esta calculado para D=64 (NTT, matrices, cotas);
 * en D=256 se ESPERA que falle. Es una sonda, no un parametro soportado. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fastrandombytes.h"
#include "gentdrow.h"
#include "keygen.h"
#include "mint.h"
#include "param.h"
#include "poly_q.h"
#include "poly_qbig.h"
#include "randombytes.h"
#include "sammat.h"
#include "spend.h"
#include "verify.h"

static unsigned char seed[CRYPTO_BYTES];
static POLY_Q g[M + R][N], h[M];
static POLY_QBIG g_hat[M_HAT + V_HAT][N_HAT];
static POLY_Q sk[N_SPENT][M], ck[N_SPENT][M], pk_out[S][N];
static SPEND_OUT out;
static POLY_Q cn_out[S][N], s;
static POLY_QBIG t0[M_HAT], t1[V_1], t2[V_2], td[N_HAT - 1];
static ACT a_in[N_SPENT];
static ASK ask;

int main(void) {
  uint64_t i, l, amt, v, amt_out[S];

  printf("=== Au ronda completa, D=%lld (el codigo original usa D=64) ===\n",
         (long long)D);
  printf("q=%lld  M=%lld  N=%lld  BETA=%lld\n", (long long)Q, (long long)M,
         (long long)N, (long long)BETA);

  srand(20261126);
  randombytes(seed, CRYPTO_BYTES);
  fastrandombytes_setseed_prv(seed);
  randombytes(seed, CRYPTO_BYTES);
  fastrandombytes_setseed_pub(seed);

  sample_mat_gr(g);
  sample_mat_gh(h);
  sample_mat_gbig(g_hat);
  gentdrow(t0, t1, t2, td, g_hat);

  l = rand() % N_SPENT;
  for (i = 0; i < N_SPENT; i++) {
    keygen(a_in[i].pk[0], sk[i], g);
    amt = (i == l) ? 1000000ULL : ((uint64_t)rand() << 32) ^ (uint64_t)rand();
    mint(a_in[i].cn[0], ck[i], amt, g);
    if (i == l) {
      ask.amt = amt;
      amt_out[0] = 456;
      amt_out[1] = amt - 456;
    }
  }
  for (i = 0; i < M; i++) {
    memcpy(ask.r + i, sk[l] + i, sizeof(POLY_Q));
    memcpy(ask.ck + i, ck[l] + i, sizeof(POLY_Q));
  }
  for (i = 0; i < N; i++) {
    memcpy(pk_out[0] + i, a_in[0].pk[0] + i, sizeof(POLY_Q));
    memcpy(pk_out[1] + i, a_in[1].pk[0] + i, sizeof(POLY_Q));
  }

  spend(cn_out, &s, &out, a_in, l, &ask, pk_out, amt_out, g, h, g_hat, t0, t1,
        t2);
  v = verify(a_in, pk_out, cn_out, &out, &s, g, h, g_hat, t0, t1, t2);

  printf("verificacion: %s (v=%llu)\n", v ? "OK" : "FALLA",
         (unsigned long long)v);
  printf("pago=%llu/%llu\n", (unsigned long long)amt_out[0],
         (unsigned long long)amt_out[1]);
  return v ? 0 : 1;
}
