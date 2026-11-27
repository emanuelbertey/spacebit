/* Au D=256 en modo BLOQUES (sin NTT densa). Firma + verifica y mide tiempo. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "fastrandombytes.h"
#include "gentdrow.h"
#include "keygen.h"
#include "mint.h"
#include "param.h"
#include "randombytes.h"
#include "sammat.h"
#include "spend.h"
#include "verify.h"

extern uint64_t qbig_mode;

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
  clock_t c0, c1, c2;

  qbig_mode = 1;  /* solo subanillos */

  printf("=== Au D=%lld, bloques, llaves=%lld ===\n", (long long)D,
         (long long)N_SPENT);
  fflush(stdout);

  srand(20261118);
  randombytes(seed, CRYPTO_BYTES);
  fastrandombytes_setseed_prv(seed);
  randombytes(seed, CRYPTO_BYTES);
  fastrandombytes_setseed_pub(seed);

  sample_mat_gr(g);
  sample_mat_gh(h);
  sample_mat_gbig(g_hat);
  gentdrow(t0, t1, t2, td, g_hat);

  l = rand() % N_SPENT;
  c0 = clock();
  for (i = 0; i < N_SPENT; i++) {
    keygen(a_in[i].pk[0], sk[i], g);
    amt = (i == l) ? 1000000ULL
                   : (((uint64_t)rand() << 32) ^ (uint64_t)rand());
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
  c1 = clock();
  v = verify(a_in, pk_out, cn_out, &out, &s, g, h, g_hat, t0, t1, t2);
  c2 = clock();

  printf("spend=%.3f s verify=%.3f s\n",
         (double)(c1 - c0) / CLOCKS_PER_SEC,
         (double)(c2 - c1) / CLOCKS_PER_SEC);
  printf("verificacion: %s (v=%llu)\n", v ? "OK" : "FALLA",
         (unsigned long long)v);
  printf("pago=%llu/%llu\n", (unsigned long long)amt_out[0],
         (unsigned long long)amt_out[1]);
  fflush(stdout);
  return v ? 0 : 1;
}
