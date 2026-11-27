/* Mide la distribucion REAL de z_c en el Au portado a D=128.
 * No cambia ningun parametro: solo observa las respuestas que produce spend().
 * z_c = x*rc + rd, cota del verificador ZC_NORM = BBIG - 448.
 * Devuelve 0 siempre; imprime estadisticas en formato JSON.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fastrandombytes.h"
#include "gentdrow.h"
#include "keygen.h"
#include "mint.h"
#include "param.h"
#include "poly_q.h"
#include "randombytes.h"
#include "sammat.h"
#include "spend.h"
#include "verify.h"

#define ROUNDS 8
#define NSAMPLES (ROUNDS * M * D)

static unsigned char semilla[CRYPTO_BYTES];
static POLY_Q g[M + R][N], h[M];
static POLY_QBIG g_hat[M_HAT + V_HAT][N_HAT];
static POLY_Q sk[N_SPENT][M], ck[N_SPENT][M], pk_out[S][N];
static SPEND_OUT out;
static POLY_Q cn_out[S][N], s;
static POLY_QBIG t0[M_HAT], t1[V_1], t2[V_2], td[N_HAT - 1];
static ACT a_in[N_SPENT];
static ASK ask;
static uint64_t muestras[NSAMPLES];

static int cmp_u64(const void *a, const void *b) {
  uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
  return (x > y) - (x < y);
}

int main(void) {
  uint64_t r, i, j, pos = 0, v, absv, max = 0;
  uint64_t sobre = 0, mitad = 0, cuarto = 0, octavo = 0;
  uint64_t p50, p90, p99, p999, p9999;
  uint64_t amt, amt_out[S], l;

  srand(20261118);
  randombytes(semilla, CRYPTO_BYTES);
  fastrandombytes_setseed_prv(semilla);
  randombytes(semilla, CRYPTO_BYTES);
  fastrandombytes_setseed_pub(semilla);

  sample_mat_gr(g);
  sample_mat_gh(h);
  sample_mat_gbig(g_hat);
  gentdrow(t0, t1, t2, td, g_hat);

  for (r = 0; r < ROUNDS; r++) {
    l = rand() % N_SPENT;
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
    spend(cn_out, &s, &out, a_in, l, &ask, pk_out, amt_out, g, h, g_hat, t0,
          t1, t2);

    for (i = 0; i < M; i++) {
      for (j = 0; j < D; j++) {
        /* zc se guarda como int64 con signo (complemento a dos), no mod Q */
        v = out.zc[i].poly[j];
        absv = (v <= (uint64_t)INT64_MAX) ? v : (uint64_t)(-(int64_t)v);
        muestras[pos++] = absv;
        if (absv > max) max = absv;
        if (absv > (uint64_t)ZC_NORM) sobre++;
        if (absv > (uint64_t)ZC_NORM / 2) mitad++;
        if (absv > (uint64_t)ZC_NORM / 4) cuarto++;
        if (absv > (uint64_t)ZC_NORM / 8) octavo++;
      }
    }
  }

  qsort(muestras, NSAMPLES, sizeof(uint64_t), cmp_u64);
  p50 = muestras[NSAMPLES * 50 / 100];
  p90 = muestras[NSAMPLES * 90 / 100];
  p99 = muestras[NSAMPLES * 99 / 100];
  p999 = muestras[NSAMPLES * 999 / 1000];
  p9999 = muestras[NSAMPLES * 9999 / 10000];

  printf("{\n");
  printf("  \"esquema\": \"au\", \"dim\": %lld, \"rounds\": %d,\n",
         (long long)D, ROUNDS);
  printf("  \"muestras\": %llu,\n", (unsigned long long)NSAMPLES);
  printf("  \"ZC_NORM\": %lld, \"log2_ZC_NORM\": %.2f,\n",
         (long long)ZC_NORM, (double)log2((double)ZC_NORM));
  printf("  \"max_abs\": %llu, \"log2_max\": %.2f,\n",
         (unsigned long long)max, (double)log2((double)(max ? max : 1)));
  printf("  \"p50\": %llu, \"p90\": %llu, \"p99\": %llu, "
         "\"p99_9\": %llu, \"p99_99\": %llu,\n",
         (unsigned long long)p50, (unsigned long long)p90,
         (unsigned long long)p99, (unsigned long long)p999,
         (unsigned long long)p9999);
  printf("  \"sobre_ZC_NORM\": %llu, \"sobre_mitad\": %llu, "
         "\"sobre_cuarto\": %llu, \"sobre_octavo\": %llu,\n",
         (unsigned long long)sobre, (unsigned long long)mitad,
         (unsigned long long)cuarto, (unsigned long long)octavo);
  printf("  \"margen_bits\": %.2f\n",
         (double)log2((double)ZC_NORM) - (double)log2((double)(max ? max : 1)));
  printf("}\n");
  return 0;
}
