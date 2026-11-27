/* Una ronda del Plus (D=256) con N_SPENT llaves (50..100 por compile -D).
 * Comprueba que firma, verifica y que los montos cuadran.
 * Sin getchar(): devuelve 0 si todo OK, 1 si falla.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "cpucycles.h"
#include "fastrandombytes.h"
#include "keygen.h"
#include "mint.h"
#include "param.h"
#include "poly_q.h"
#include "randombytes.h"
#include "sammat.h"
#include "setup.h"
#include "spend.h"
#include "verify.h"

static unsigned char seed[CRYPTO_BYTES];
static POLY_Q mat1[K + 3][N], mat2[K][3];
static POLY_QHAT mat3[K_HAT + 2 * N_SPENT + 4][N_HAT];

int main(void) {
  static ACT a_in[N_SPENT][M];
  static ASK ask[M];
  static PK_OUT pk_out[S];
  static uint64_t amt_out[S];
  static POLY_R sk[N_SPENT][M][N_BAR + K];
  static unsigned char sn[N_SPENT][M][CRYPTO_BYTES];
  static POLY_R ck[N_SPENT][M][N + K + 1];
  static uint64_t a_hat[N_SPENT][M][R];
  static SPEND_OUT spend_out;
  static CN_OUT cn_out[S];
  static CK_OUT ck_out[S];
  static SN_OUT sn_out[M];
  static unsigned char pp[CRYPTO_BYTES];

  uint64_t i, l, amt0, amt1, v, total_in;
  clock_t c0, c1;

  printf("=== PLUS D=256, llaves=%lld, BG_0=%lld ===\n",
         (long long)N_SPENT, (long long)BG_0);
  fflush(stdout);

  srand(20261118);
  randombytes(seed, CRYPTO_BYTES);
  fastrandombytes_setseed_prv(seed);
  randombytes(seed, CRYPTO_BYTES);
  fastrandombytes_setseed_pub(seed);
  setup(pp, seed);

  sample_mat1(mat1);
  sample_mat2(mat2);
  sample_mat3(mat3);

  l = rand() % N_SPENT;
  c0 = clock();
  for (i = 0; i < N_SPENT; i++) {
    keygen(a_in[i][0].pk, sk[i][0], sn[i][0], mat1);
    keygen(a_in[i][1].pk, sk[i][1], sn[i][1], mat1);
    amt0 = (((uint64_t)rand() << 32) + rand()) % (1ULL << 40);
    amt1 = (((uint64_t)rand() << 32) + rand()) % (1ULL << 40);
    mint(a_in[i][0].cn, ck[i][0], a_hat[i][0], amt0, mat1, mat2);
    mint(a_in[i][1].cn, ck[i][1], a_hat[i][1], amt1, mat1, mat2);
    if (i == l) {
      ask[0].amt_in = amt0;
      ask[1].amt_in = amt1;
      amt_out[0] = rand() % (amt0 + amt1);
      amt_out[1] = amt0 + amt1 - amt_out[0];
      total_in = amt0 + amt1;
    }
  }
  for (i = 0; i < N_BAR + K; i++) {
    memcpy(ask[0].sk + i, sk[l][0] + i, sizeof(POLY_R));
    memcpy(ask[1].sk + i, sk[l][1] + i, sizeof(POLY_R));
  }
  memcpy(ask[0].sn, sn[l][0], CRYPTO_BYTES);
  memcpy(ask[1].sn, sn[l][1], CRYPTO_BYTES);
  for (i = 0; i < N + K + 1; i++) {
    memcpy(ask[0].ck + i, ck[l][0] + i, sizeof(POLY_R));
    memcpy(ask[1].ck + i, ck[l][1] + i, sizeof(POLY_R));
  }
  for (i = 0; i < N_BAR; i++) {
    memcpy(pk_out[0].pk + i, a_in[0][0].pk + i, sizeof(POLY_Q));
    memcpy(pk_out[1].pk + i, a_in[0][1].pk + i, sizeof(POLY_Q));
  }
  spend(&spend_out, cn_out, ck_out, sn_out, a_in, l, ask, pk_out, amt_out,
        mat1, mat2, mat3, pp);
  c1 = clock();

  v = verify(a_in, pk_out, cn_out, sn_out, &spend_out, mat1, mat2, mat3, pp);

  printf("firmante=%llu  entrada=%llu  pago=%llu  vuelto=%llu\n",
         (unsigned long long)l, (unsigned long long)total_in,
         (unsigned long long)amt_out[0], (unsigned long long)amt_out[1]);
  printf("montos cuadran: %s\n",
         (amt_out[0] + amt_out[1] == total_in) ? "OK" : "FALLA");
  printf("tiempo firmar: %.3f s\n",
         (double)(c1 - c0) / (double)CLOCKS_PER_SEC);
  printf("verificacion: %s (v=%llu)\n", v ? "OK" : "FALLA",
         (unsigned long long)v);
  fflush(stdout);

  if (!v) return 1;
  if (amt_out[0] + amt_out[1] != total_in) return 1;
  return 0;
}
