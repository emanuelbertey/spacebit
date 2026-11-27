/* Au D=256 en modo BLOQUES (sin NTT densa). Firma + verifica y mide tiempo. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "audit.h"
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

#define LOG(...) do { printf(__VA_ARGS__); printf("\n"); fflush(stdout); } while (0)

int main(void) {
  uint64_t i, l, amt, v, amt_out[S], l_au, amt_out_au[S], aud_ok;
  clock_t c0, c1, c2;

  qbig_mode = 1;  /* solo subanillos */

  LOG("=== Au D=%lld, bloques, llaves=%lld | M=%lld rank=%lld ===",
      (long long)D, (long long)N_SPENT, (long long)M,
      (long long)(M / 19));

  srand(20261118);
  randombytes(seed, CRYPTO_BYTES);
  fastrandombytes_setseed_prv(seed);
  randombytes(seed, CRYPTO_BYTES);
  fastrandombytes_setseed_pub(seed);

  LOG("[1/6] sample_mat...");
  sample_mat_gr(g);
  sample_mat_gh(h);
  sample_mat_gbig(g_hat);
  LOG("[2/6] gentdrow...");
  gentdrow(t0, t1, t2, td, g_hat);

  l = rand() % N_SPENT;
  LOG("[3/6] keygen+mint (%d llaves, firmante=%d)...", (int)N_SPENT, (int)l);
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
  LOG("[4/6] spend (rechazo: ver [spend] en el log)...");
  spend(cn_out, &s, &out, a_in, l, &ask, pk_out, amt_out, g, h, g_hat, t0, t1,
        t2);
  c1 = clock();
  LOG("[5/6] verify...");
  v = verify(a_in, pk_out, cn_out, &out, &s, g, h, g_hat, t0, t1, t2);
  c2 = clock();

  LOG("spend=%.3f s verify=%.3f s",
      (double)(c1 - c0) / CLOCKS_PER_SEC,
      (double)(c2 - c1) / CLOCKS_PER_SEC);
  LOG("verificacion: %s (v=%llu)", v ? "OK" : "FALLA", (unsigned long long)v);

  /* auditoria Au: recupera el firmante y los montos desde g_hat y td */
  LOG("[6/6] audit (llave/montos)...");
  audit(&l_au, amt_out_au, &out, g_hat, td);
  aud_ok = (l_au == l) && (amt_out_au[0] == amt_out[0]) &&
           (amt_out_au[1] == amt_out[1]);
  LOG("auditoria (llave/montos): %s (l_au=%llu, pago_au=%llu/%llu)",
      aud_ok ? "OK" : "FALLA", (unsigned long long)l_au,
      (unsigned long long)amt_out_au[0], (unsigned long long)amt_out_au[1]);
  LOG("pago=%llu/%llu", (unsigned long long)amt_out[0],
      (unsigned long long)amt_out[1]);
  return (v && aud_ok) ? 0 : 1;
}
