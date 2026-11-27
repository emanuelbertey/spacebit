/* Prueba REAL de Au con D=128 (en vez de 64).
 * Ronda completa: keygen + mint + spend + verify, mas autotests de NTT.
 * NO toca fuentes originales: se compila sobre `auditor_d128/` generado por
 * test-d128/port.py en GitHub Actions.
 * Devuelve 0 si todo verifica, 1 si la ronda falla, 2 si falla un NTT.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fastrandombytes.h"
#include "gentdrow.h"
#include "keygen.h"
#include "mint.h"
#include "param.h"
#include "poly_mult.h"
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

/* Comprueba que el NTT de qbig (reemplazado en el port) multiplica bien:
 * compara contra la convolucion negaciclica escolar. */
static int selftest_qbig(void) {
  static POLY_QBIG a, b, c, d, ref;
  uint64_t i, j, k, p1, p2;

  for (i = 0; i < D; i++) {
    a.poly1[i] = (i * 1000003ULL + 12345ULL) % QBIG1;
    a.poly2[i] = (i * 2000003ULL + 54321ULL) % QBIG2;
    b.poly1[i] = (i * 3000017ULL + 777ULL) % QBIG1;
    b.poly2[i] = (i * 4000037ULL + 999ULL) % QBIG2;
    ref.poly1[i] = 0;
    ref.poly2[i] = 0;
  }

  for (i = 0; i < D; i++) {
    for (j = 0; j < D; j++) {
      p1 = (a.poly1[i] * b.poly1[j]) % QBIG1;
      p2 = (a.poly2[i] * b.poly2[j]) % QBIG2;
      k = i + j;
      if (k < D) {
        ref.poly1[k] = (ref.poly1[k] + p1) % QBIG1;
        ref.poly2[k] = (ref.poly2[k] + p2) % QBIG2;
      } else {
        ref.poly1[k - D] = (ref.poly1[k - D] + QBIG1 - p1) % QBIG1;
        ref.poly2[k - D] = (ref.poly2[k - D] + QBIG2 - p2) % QBIG2;
      }
    }
  }

  c = a;
  d = b;
  ntt_qbig(&c);
  ntt_qbig(&d);
  mult_rqbig(&c, &c, &d);
  intt_qbig(&c);

  for (i = 0; i < D; i++) {
    if (c.poly1[i] != ref.poly1[i] || c.poly2[i] != ref.poly2[i]) {
      printf("NTT qbig selftest: FALLA en %llu (%llu/%llu vs %llu/%llu)\n",
             (unsigned long long)i, (unsigned long long)c.poly1[i],
             (unsigned long long)c.poly2[i], (unsigned long long)ref.poly1[i],
             (unsigned long long)ref.poly2[i]);
      return 1;
    }
  }
  printf("NTT qbig selftest: OK\n");
  return 0;
}

/* Lo mismo para R_q (subanillos de grado D/4). */
static int selftest_rq(void) {
  static POLY_Q a, b, c, ref;
  uint64_t i, j, k, p;

  for (i = 0; i < D; i++) {
    a.poly[i] = (i * 7 + 3) % 1000;
    b.poly[i] = (i * 13 + 5) % 1000;
    ref.poly[i] = 0;
  }

  for (i = 0; i < D; i++) {
    for (j = 0; j < D; j++) {
      p = (a.poly[i] * b.poly[j]) % Q;
      k = i + j;
      if (k < D) {
        ref.poly[k] = (ref.poly[k] + p) % Q;
      } else {
        ref.poly[k - D] = (ref.poly[k - D] + Q - p) % Q;
      }
    }
  }

  ntt_q(&a);
  ntt_q(&b);
  mult_rq(&c, &a, &b);
  intt_q(&c);

  for (i = 0; i < D; i++) {
    if (c.poly[i] != ref.poly[i]) {
      printf("NTT rq selftest: FALLA en %llu (%llu vs %llu)\n",
             (unsigned long long)i, (unsigned long long)c.poly[i],
             (unsigned long long)ref.poly[i]);
      return 1;
    }
  }
  printf("NTT rq selftest: OK\n");
  return 0;
}

int main(void) {
  uint64_t i, l, amt, v, amt_out[S];

  printf("=== Au ronda completa, D=%lld (original D=64) ===\n", (long long)D);
  printf("q=%lld  M=%lld  N=%lld  BETA=%lld  W=%lld\n", (long long)Q,
         (long long)M, (long long)N, (long long)BETA, (long long)W);
  fflush(stdout);

  if (selftest_qbig() || selftest_rq())
    return 2;

  srand(20261126);
  randombytes(seed, CRYPTO_BYTES);
  fastrandombytes_setseed_prv(seed);
  randombytes(seed, CRYPTO_BYTES);
  fastrandombytes_setseed_pub(seed);

  sample_mat_gr(g);
  sample_mat_gh(h);
  sample_mat_gbig(g_hat);
  gentdrow(t0, t1, t2, td, g_hat);
  printf("setup (g, h, g_hat, gadget): OK\n");
  fflush(stdout);

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
  printf("keygen+mint de %lld llaves: OK (gasto en %llu)\n",
         (long long)N_SPENT, (unsigned long long)l);
  fflush(stdout);

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
  printf("spend: OK\n");
  fflush(stdout);

  v = verify(a_in, pk_out, cn_out, &out, &s, g, h, g_hat, t0, t1, t2);

  printf("verificacion: %s (v=%llu)\n", v ? "OK" : "FALLA",
         (unsigned long long)v);
  printf("pago=%llu/%llu\n", (unsigned long long)amt_out[0],
         (unsigned long long)amt_out[1]);
  fflush(stdout);
  return v ? 0 : 1;
}
