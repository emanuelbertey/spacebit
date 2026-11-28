#!/usr/bin/env python3
"""Porta Au de D=64 a D=128 generando una copia `auditor_d128/`.

D=128 es el maximo que permiten los modulos actuales: el NTT de qbig necesita
una raiz primitiva de orden 2D, y QBIG2-1 solo tiene 2^8, asi que 2D <= 256
=> D <= 128. D=256 exigiria cambiar QBIG2 (y rehacer la seguridad).

No toca `auditor/` (original intacto). Recalcula:
  - D=128 en poly_param.h
  - todos los largos de muestreo (REJ_*_LEN, SAMPLE_*_LEN) y HASH_OUTPUT_LEN
  - la mascara de parse_hash (& 63 -> & (D-1))
  - la multiplicacion schoolbook de subanillo (antes fija a grado 16) a D/4=32
  - el NTT de qbig: se reemplaza poly_qbig.c por una implementacion generica
    (aritmetica modular simple, raiz de orden 2D calculada en runtime)
  - las tablas del gadget (gentdrow.c) a partir de los polinomios recuperados
    con la implementacion ORIGINAL de D=64 (archivo --gadget)
  - TG opcionalmente escalado (--tg-factor); por defecto 1 para que la
    auditoria estatica siga viendo los valores de fabrica.

Todo se ejecuta en GitHub Actions; este script no compila nada.
"""
import argparse
import pathlib
import re
import shutil
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
SRC = ROOT / "auditor"
DST = ROOT / "auditor_d128"

D = 128
D_OLD = 64
Q = 2147221513
QBIG1 = 134215681
QBIG2 = 536870657

# Parametros que NO dependen de D.
N, M, R = 18, 38, 64
N_HAT, M_HAT, BETA, S, W = 35, 69, 100, 2, 56
V_1 = V_2 = BETA + S * R + R - 1          # 291
V_HAT = 2 * (BETA + R - 1 + S * R)        # 582
TG_NAMES = ("TG_COEFF_0", "TG_COEFF_1", "TG_COEFF_2", "TG_COEFF_3")


def slen(n):
    """Muestras necesarias + margen (sobre-asignar es seguro)."""
    return int(n * 1.15) + 64


def setdef(text, name, value):
    pat = re.compile(r"(#define\s+" + re.escape(name) + r"\s+)\S+")
    new, cnt = pat.subn(lambda m: m.group(1) + str(value), text, count=1)
    if cnt != 1:
        raise SystemExit(f"no pude reemplazar #define {name}")
    return new


def getdef(text, name):
    m = re.search(r"#define\s+" + re.escape(name) + r"\s+(\S+)", text)
    if not m:
        raise SystemExit(f"no encontre #define {name}")
    return int(m.group(1).replace("LL", "").replace("U", ""))


# ---------------------------------------------------------------- NTT limpio
def modpow(b, e, q):
    return pow(b, e, q)


def root_2d(q, d):
    """Primera raiz de orden exacto 2d modulo q (mismo criterio que el C)."""
    order = 2 * d
    if (q - 1) % order:
        raise SystemExit(f"orden {order} no divide q-1 para q={q}")
    e = (q - 1) // order
    g = 2
    while g < 1000000:
        r = pow(g, e, q)
        if pow(r, order, q) == 1 and pow(r, order // 2, q) != 1:
            return r
        g += 1
    raise SystemExit(f"no encontre raiz de orden {order} mod {q}")


def _bitrev(a):
    n = len(a)
    j = 0
    for i in range(1, n):
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j ^= bit
        if i < j:
            a[i], a[j] = a[j], a[i]


def ntt_py(a, q, psi):
    """NTT negaciclica natural: A[k]=sum_j a[j]*psi^((2k+1)j)."""
    n = len(a)
    w = psi * psi % q
    pw = 1
    for j in range(n):
        a[j] = (a[j] % q) * pw % q
        pw = pw * psi % q
    _bitrev(a)
    length = 2
    while length <= n:
        wl = pow(w, n // length, q)
        for i in range(0, n, length):
            ww = 1
            for j in range(length // 2):
                u = a[i + j]
                v = a[i + j + length // 2] * ww % q
                a[i + j] = (u + v) % q
                a[i + j + length // 2] = (u - v) % q
                ww = ww * wl % q
        length <<= 1
    return a


# ------------------------------------------------------------- C de poly_qbig
POLY_QBIG_C = r"""/* ****************************** *
 * ntt on R_{\hat{q}}             *
 * ****************************** *
 * Portatil para cualquier D (D=128). Aritmetica modular simple (los modulos
 * QBIG1/QBIG2 son < 2^30, asi que a*b < 2^60 entra en uint64). La raiz
 * primitiva de orden 2D se calcula en runtime y se cachea.
 */
#include "poly_qbig.h"
#include "poly_param.h"
#include <stdint.h>

static uint64_t modpow_q(uint64_t b, uint64_t e, uint64_t q) {
  uint64_t r = 1 % q;
  b %= q;
  while (e) {
    if (e & 1)
      r = (r * b) % q;
    b = (b * b) % q;
    e >>= 1;
  }
  return r;
}

static uint64_t root_2d(uint64_t q) {
  const uint64_t order = 2 * (uint64_t)D;
  const uint64_t e = (q - 1) / order;
  uint64_t g, r;
  for (g = 2; g < 1000000; g++) {
    r = modpow_q(g, e, q);
    if (modpow_q(r, order, q) == 1 && modpow_q(r, order >> 1, q) != 1)
      return r;
  }
  return 0;
}

static uint64_t psi1 = 0, psi2 = 0;

static void init_roots(void) {
  if (!psi1)
    psi1 = root_2d(QBIG1);
  if (!psi2)
    psi2 = root_2d(QBIG2);
}

static void ntt_core(uint64_t *a, const uint64_t q, const uint64_t psi) {
  const uint64_t n = (uint64_t)D;
  const uint64_t w = (psi * psi) % q;
  uint64_t i, j, len, wl, ww, u, v, pw;

  /* twist: a_j <- a_j * psi^j */
  pw = 1;
  for (j = 0; j < n; j++) {
    a[j] = (a[j] % q) * pw % q;
    pw = (pw * psi) % q;
  }

  /* bit reversal */
  for (i = 1, j = 0; i < n; i++) {
    uint64_t bit = n >> 1;
    while (j & bit) {
      j ^= bit;
      bit >>= 1;
    }
    j ^= bit;
    if (i < j) {
      u = a[i];
      a[i] = a[j];
      a[j] = u;
    }
  }

  /* Cooley-Tukey */
  for (len = 2; len <= n; len <<= 1) {
    wl = modpow_q(w, n / len, q);
    for (i = 0; i < n; i += len) {
      ww = 1;
      for (j = 0; j < (len >> 1); j++) {
        u = a[i + j];
        v = (a[i + j + (len >> 1)] * ww) % q;
        a[i + j] = (u + v) % q;
        a[i + j + (len >> 1)] = (u + q - v) % q;
        ww = (ww * wl) % q;
      }
    }
  }
}

static void intt_core(uint64_t *a, const uint64_t q, const uint64_t psi) {
  const uint64_t n = (uint64_t)D;
  const uint64_t w = (psi * psi) % q;
  const uint64_t winv = modpow_q(w, q - 2, q);
  const uint64_t ninv = modpow_q(n, q - 2, q);
  const uint64_t psi_inv = modpow_q(psi, q - 2, q);
  uint64_t i, j, len, wl, ww, u, v, pw;

  for (j = 0; j < n; j++)
    a[j] %= q;

  /* bit reversal */
  for (i = 1, j = 0; i < n; i++) {
    uint64_t bit = n >> 1;
    while (j & bit) {
      j ^= bit;
      bit >>= 1;
    }
    j ^= bit;
    if (i < j) {
      u = a[i];
      a[i] = a[j];
      a[j] = u;
    }
  }

  /* Cooley-Tukey con winv: MISMA mariposa del forward (ww multiplica a v),
   * no la mariposa de Gentleman-Sande. Mezclarlas rompe el roundtrip. */
  for (len = 2; len <= n; len <<= 1) {
    wl = modpow_q(winv, n / len, q);
    for (i = 0; i < n; i += len) {
      ww = 1;
      for (j = 0; j < (len >> 1); j++) {
        u = a[i + j];
        v = (a[i + j + (len >> 1)] * ww) % q;
        a[i + j] = (u + v) % q;
        a[i + j + (len >> 1)] = (u + q - v) % q;
        ww = (ww * wl) % q;
      }
    }
  }

  /* escala por n^{-1} y deshace el twist */
  pw = 1;
  for (j = 0; j < n; j++) {
    a[j] = a[j] * ninv % q;
    a[j] = a[j] * pw % q;
    pw = (pw * psi_inv) % q;
  }
}

void ntt_qbig(POLY_QBIG *a) {
  uint64_t i;
  init_roots();
  for (i = 0; i < D; i++) {
    a->poly1[i] %= QBIG1;
    a->poly2[i] %= QBIG2;
  }
  ntt_core(a->poly1, QBIG1, psi1);
  ntt_core(a->poly2, QBIG2, psi2);
}

void intt_qbig(POLY_QBIG *a) {
  init_roots();
  intt_core(a->poly1, QBIG1, psi1);
  intt_core(a->poly2, QBIG2, psi2);
}

void mult_rqbig(POLY_QBIG *out, const POLY_QBIG *a, const POLY_QBIG *b) {
  uint64_t i;
  for (i = 0; i < D; i++) {
    out->poly1[i] = (a->poly1[i] * b->poly1[i]) % QBIG1;
    out->poly2[i] = (a->poly2[i] * b->poly2[i]) % QBIG2;
  }
}
"""


# Implementacion alternativa estilo Plus: D = K*S con S=4. Solo necesita
# raices de orden 2K (128 para D=256), no 2D (512). La representacion es
# a[r*K+i] = A_r(eta^(2i+1)), con A_r(Y)=sum_q a[q*S+r]Y^q.
POLY_QBIG_SUB_C = r"""/* Multiplicacion en qbig por subanillos (estilo Plus).
 * D = K*S, S=4. ntt/intt son DFT de longitud K sobre cada residuo r mod S;
 * mult es schoolbook de grado S con vuelta multiplicada por eta^(2i+1).
 * Raiz requerida: orden 2K = 2D/S (128 con D=256), no 2D (512).
 */
#include "poly_qbig.h"
#include "poly_param.h"
#include <stdint.h>

#define SUB_S 4
#define SUB_K (D / SUB_S)

uint64_t qbig_mode = 0; /* 0 = NTT densa, 1 = subanillos */

extern void ntt_qbig_dense(POLY_QBIG *a);
extern void intt_qbig_dense(POLY_QBIG *a);
extern void mult_rqbig_dense(POLY_QBIG *out, const POLY_QBIG *a,
                             const POLY_QBIG *b);

static uint64_t mpow_s(uint64_t b, uint64_t e, uint64_t q) {
  uint64_t r = 1 % q;
  b %= q;
  while (e) {
    if (e & 1) r = r * b % q;
    b = b * b % q;
    e >>= 1;
  }
  return r;
}

static uint64_t root_order_s(uint64_t q, uint64_t order) {
  uint64_t e = (q - 1) / order, g, r;
  for (g = 2; g < 1000000; g++) {
    r = mpow_s(g, e, q);
    if (mpow_s(r, order, q) == 1 && mpow_s(r, order >> 1, q) != 1)
      return r;
  }
  return 0;
}

static uint64_t eta1 = 0, eta2 = 0;

static void init_sub(void) {
  if (!eta1) eta1 = root_order_s(QBIG1, 2 * (uint64_t)SUB_K);
  if (!eta2) eta2 = root_order_s(QBIG2, 2 * (uint64_t)SUB_K);
}

static void bitrev_s(uint64_t *a, uint64_t n) {
  uint64_t i, j = 0, bit, u;
  for (i = 1; i < n; i++) {
    bit = n >> 1;
    while (j & bit) { j ^= bit; bit >>= 1; }
    j ^= bit;
    if (i < j) { u = a[i]; a[i] = a[j]; a[j] = u; }
  }
}

static void ct_core_s(uint64_t *a, uint64_t q, uint64_t w, uint64_t n) {
  uint64_t i, j, len, wl, ww, u, v;
  bitrev_s(a, n);
  for (len = 2; len <= n; len <<= 1) {
    wl = mpow_s(w, n / len, q);
    for (i = 0; i < n; i += len) {
      ww = 1;
      for (j = 0; j < (len >> 1); j++) {
        u = a[i + j];
        v = a[i + j + (len >> 1)] * ww % q;
        a[i + j] = (u + v) % q;
        a[i + j + (len >> 1)] = (u + q - v) % q;
        ww = ww * wl % q;
      }
    }
  }
}

/* Ojo: no se puede escribir en `a` mientras se siguen leyendo los demas
 * residuos r (los indices qi*SUB_S+r se pisan). Se guarda todo en `buf`. */
static void ntt_sub_core(uint64_t *a, uint64_t q, uint64_t eta) {
  uint64_t buf[SUB_S][SUB_K];
  uint64_t zeta = eta * eta % q;
  uint64_t r, qi, i, pw;
  for (r = 0; r < SUB_S; r++) {
    pw = 1;
    for (qi = 0; qi < SUB_K; qi++) {
      buf[r][qi] = (a[qi * SUB_S + r] % q) * pw % q;
      pw = pw * eta % q;
    }
    ct_core_s(buf[r], q, zeta, SUB_K);
  }
  for (r = 0; r < SUB_S; r++)
    for (i = 0; i < SUB_K; i++) a[r * SUB_K + i] = buf[r][i];
}

static void intt_sub_core(uint64_t *a, uint64_t q, uint64_t eta) {
  uint64_t buf[SUB_S][SUB_K];
  uint64_t zeta_inv = mpow_s(eta * eta % q, q - 2, q);
  uint64_t eta_inv = mpow_s(eta, q - 2, q);
  uint64_t kinv = mpow_s(SUB_K, q - 2, q);
  uint64_t r, qi, i, pw;
  for (r = 0; r < SUB_S; r++) {
    for (i = 0; i < SUB_K; i++) buf[r][i] = a[r * SUB_K + i] % q;
    ct_core_s(buf[r], q, zeta_inv, SUB_K);
    pw = 1;
    for (qi = 0; qi < SUB_K; qi++) {
      buf[r][qi] = buf[r][qi] * kinv % q * pw % q;
      pw = pw * eta_inv % q;
    }
  }
  for (r = 0; r < SUB_S; r++)
    for (qi = 0; qi < SUB_K; qi++) a[qi * SUB_S + r] = buf[r][qi];
}

static void mult_sub_core(uint64_t *out, const uint64_t *a, const uint64_t *b,
                          uint64_t q, uint64_t eta) {
  uint64_t i, j, k, r, beta, prod;
  uint64_t ca[SUB_S], cb[SUB_S], co[SUB_S];
  for (i = 0; i < SUB_K; i++) {
    beta = mpow_s(eta, 2 * i + 1, q);
    for (r = 0; r < SUB_S; r++) {
      ca[r] = a[r * SUB_K + i] % q;
      cb[r] = b[r * SUB_K + i] % q;
      co[r] = 0;
    }
    for (j = 0; j < SUB_S; j++) {
      for (k = 0; k < SUB_S; k++) {
        prod = ca[j] * cb[k] % q;
        if (j + k < SUB_S)
          co[j + k] = (co[j + k] + prod) % q;
        else
          co[j + k - SUB_S] = (co[j + k - SUB_S] + beta * prod) % q;
      }
    }
    for (r = 0; r < SUB_S; r++) out[r * SUB_K + i] = co[r];
  }
}

void ntt_qbig_sub(POLY_QBIG *a) {
  init_sub();
  ntt_sub_core(a->poly1, QBIG1, eta1);
  ntt_sub_core(a->poly2, QBIG2, eta2);
}

void intt_qbig_sub(POLY_QBIG *a) {
  init_sub();
  intt_sub_core(a->poly1, QBIG1, eta1);
  intt_sub_core(a->poly2, QBIG2, eta2);
}

void mult_rqbig_sub(POLY_QBIG *out, const POLY_QBIG *a, const POLY_QBIG *b) {
  init_sub();
  mult_sub_core(out->poly1, a->poly1, b->poly1, QBIG1, eta1);
  mult_sub_core(out->poly2, a->poly2, b->poly2, QBIG2, eta2);
}

void ntt_qbig(POLY_QBIG *a) {
  if (qbig_mode) ntt_qbig_sub(a); else ntt_qbig_dense(a);
}

void intt_qbig(POLY_QBIG *a) {
  if (qbig_mode) intt_qbig_sub(a); else intt_qbig_dense(a);
}

void mult_rqbig(POLY_QBIG *out, const POLY_QBIG *a, const POLY_QBIG *b) {
  if (qbig_mode) mult_rqbig_sub(out, a, b);
  else mult_rqbig_dense(out, a, b);
}
"""


MULT_OLD_RE = re.compile(
    r"/\* out=a\*b on 4 ntt subrings \*/\s*"
    r"static inline void mult\(uint64_t \*out, const uint64_t \*a, "
    r"const uint64_t \*b\)\s*\{.*?\n\}",
    re.S,
)

MULT_NEW = r"""/* out=a*b on 4 ntt subrings, grado D/4 (generico para D=128) */
static inline void mult_sub_g(uint64_t *out, const uint64_t *a,
                              const uint64_t *b, const uint64_t omega,
                              const int plus) {
  const uint64_t K = D >> 2;
  uint64_t i, j, acc, wr;
  for (i = 0; i < K; i++) {
    acc = 0;
    for (j = 0; j <= i; j++)
      acc += montgomery_q(a[j], b[i - j]);
    wr = 0;
    for (j = i + 1; j < K; j++)
      wr += montgomery_q(a[j], b[K + i - j]);
    wr = montgomery_q(red_short_q(wr), omega);
    out[i] = red_short_q(plus ? (acc + wr) : (acc + Q - wr));
  }
}

static inline void mult(uint64_t *out, const uint64_t *a, const uint64_t *b) {
  mult_sub_g(out, a, b, omega[1], 1);
  mult_sub_g(out + (D >> 2), a + (D >> 2), b + (D >> 2), omega[1], 0);
  mult_sub_g(out + (D >> 1), a + (D >> 1), b + (D >> 1), omega[3], 1);
  mult_sub_g(out + (D >> 1) + (D >> 2), a + (D >> 1) + (D >> 2),
             b + (D >> 1) + (D >> 2), omega[3], 0);
}"""


# --------------------------------------------------------- tablas de gentdrow
def c_poly_entry(p1, p2):
    def fmt(vals):
        rows = []
        for i in range(0, len(vals), 8):
            rows.append("      " + ", ".join(str(v) for v in vals[i:i + 8]))
        return ",\n".join(rows)

    return "    {{%s},\n     {%s}}" % (fmt(p1), fmt(p2))


def table_block(decl, entries):
    body = ",\n".join(c_poly_entry(p1, p2) for p1, p2 in entries)
    return decl + "\n" + body + "\n};\n"


def replace_table(text, decl_start, new_block):
    # Ojo: la tabla original cierra con "}};" en la misma linea de la ultima
    # entrada, no con un "};" suelto. Delimitamos por la siguiente declaracion.
    lines = text.splitlines()
    start = None
    for i, line in enumerate(lines):
        if line.startswith(decl_start):
            start = i
            break
    if start is None:
        raise SystemExit(f"no encontre la tabla {decl_start}")
    end = None
    for i in range(start + 1, len(lines)):
        if (lines[i].startswith("static const POLY_QBIG tg_") or
                lines[i].startswith("void gentdrow(")):
            end = i
            break
    if end is None:
        end = len(lines)
    return "\n".join(lines[:start] + new_block.splitlines() + lines[end:]) + "\n"


T_BAR_QBIG1 = 113567114
T_BAR_QBIG2 = 82595485
R1 = (1 << 32) % QBIG1
R2 = (1 << 32) % QBIG2


def load_gadget(path):
    raw = pathlib.Path(path).read_text(encoding="utf-8").split()
    vals = [int(x) for x in raw]
    need = (N_SPENT_TABLES + S * R) * D_OLD * 2
    print(f"gadget: {len(vals)} enteros (esperaba {need})")
    vals += [0] * max(0, need - len(vals))

    entries = []
    pos = 0
    for _ in range(N_SPENT_TABLES + S * R):
        p1, p2 = [], []
        for _j in range(D_OLD):
            p1.append(vals[pos] if pos < len(vals) else 0)
            p2.append(vals[pos + 1] if pos + 1 < len(vals) else 0)
            pos += 2
        entries.append((p1, p2))

    # tg_0[0] es el gadget del indice 1: t_bar*1 = [t_bar,0,...,0].
    # No asumimos si intt_qbig devuelve t_bar, t_bar*R o t_bar*R^{-1}:
    # el propio coeficiente [0] nos da la escala. Lo dividimos para obtener
    # los polinomios de bits (coeficientes 0/1) y luego los reescalamos a
    # t_bar "plano", que es lo que espera el NTT nuevo (aritmetica simple).
    s1, s2 = entries[0][0][0], entries[0][1][0]
    if s1 and s2:
        i1, i2 = pow(s1, -1, QBIG1), pow(s2, -1, QBIG2)
        bits = [([(v * i1) % QBIG1 for v in p1],
                 [(v * i2) % QBIG2 for v in p2])
                for p1, p2 in entries]
        malos = [k for k in range(min(N_SPENT_TABLES, len(bits)))
                 if any(v not in (0, 1) for v in bits[k][0])]
        if malos:
            print(f"AVISO gadget: tg_0 con valores no 0/1 en {len(malos)} "
                  f"entradas (primera {malos[0]}); se usa tal cual")
        entries = [([(v * T_BAR_QBIG1) % QBIG1 for v in p1],
                    [(v * T_BAR_QBIG2) % QBIG2 for v in p2])
                   for p1, p2 in bits]
        print(f"gadget: escala calibrada ({s1}, {s2}) -> t_bar plano")
    else:
        print("AVISO gadget: escala cero; se usan los valores recuperados tal cual")
    return entries


N_SPENT_TABLES = BETA - 1  # tg_0 tiene N_SPENT-1 entradas


def instrument_spend(text):
    """Log de intentos de rechazo en la copia generada (no toca auditor/)."""
    if "#include <stdio.h>" not in text:
        text = text.replace("#include <stdint.h>",
                            "#include <stdint.h>\n#include <stdio.h>", 1)
    decl = "uint64_t i, j, k, rej, e_norm;"
    if decl not in text:
        raise SystemExit("no encontre la declaracion de rej en spend.c")
    text = text.replace(decl, decl + "\n  static unsigned long __intento = 0;", 1)
    old = "  do {\n    rej = 0;\n"
    new = """  do {
    __intento++;
    if (__intento <= 3 || (__intento % 10) == 0)
      fprintf(stderr, "[spend] intento %lu\\n", __intento);
    if (__intento > 60) {
      fprintf(stderr, "[spend] ABORTE: 60 intentos, rechazo permanente\\n");
      exit(3);
    }
    rej = 0;
"""
    if old not in text:
        raise SystemExit("no encontre el bucle de rechazo en spend.c")
    text = text.replace(old, new, 1)

    etiquetas = ("BAP", "F00_NORM", "BRP_1", "BRP_2", "gnorm_TG",
                 "ZB_NORM", "ZC_NORM", "ZC_out0", "ZC_out1", "Z0_NORM",
                 "Z0_NORM_2", "ZM_NORM")
    viejo = "    if (rej) {\n      continue;\n    }"
    for etiqueta in etiquetas:
        if viejo not in text:
            raise SystemExit("no encontre todos los chequeos de rechazo")
        nuevo = ("    if (rej) {\n"
                 '      fprintf(stderr, "[spend] rechazo ' + etiqueta +
                 ' (intento %lu)\\n", __intento);\n'
                 "      continue;\n"
                 "    }")
        text = text.replace(viejo, nuevo, 1)
    return text


def ntt_py_sub(a, q, eta):
    """NTT por bloques (K=D/4) igual que poly_qbig_sub.c.
    Necesario para D=256, donde no existe raiz de orden 2D en QBIG2."""
    K = D >> 2
    zeta = eta * eta % q
    out = [0] * D
    for r in range(4):
        pw = 1
        tmp = []
        for qi in range(K):
            tmp.append(a[qi * 4 + r] % q * pw % q)
            pw = pw * eta % q
        _bitrev(tmp)
        length = 2
        while length <= K:
            wl = pow(zeta, K // length, q)
            for i in range(0, K, length):
                ww = 1
                for j in range(length // 2):
                    u = tmp[i + j]
                    v = tmp[i + j + length // 2] * ww % q
                    tmp[i + j] = (u + v) % q
                    tmp[i + j + length // 2] = (u - v) % q
                    ww = ww * wl % q
            length <<= 1
        for i in range(K):
            out[r * K + i] = tmp[i]
    return out


def main():
    global D, M
    ap = argparse.ArgumentParser()
    ap.add_argument("--dim", type=int, default=128)
    ap.add_argument("--n-spent", type=int, default=100)
    ap.add_argument("--M", type=int, default=None)
    ap.add_argument("--out", default="auditor_d128")
    ap.add_argument("--gadget", default=None,
                    help="gadget_polys.txt recuperado con el codigo original D=64")
    ap.add_argument("--tg-factor", type=int, default=1)
    ap.add_argument("--br-factor", type=int, default=1)
    args = ap.parse_args()
    D = args.dim
    if args.M is not None:
        M = args.M

    dst = ROOT / args.out
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(SRC, dst)

    # Logs del bucle de rechazo (solo en la copia generada).
    p = dst / "spend.c"
    p.write_text(instrument_spend(p.read_text(encoding="utf-8")), encoding="utf-8")

    # Log de tg vs norm dentro de check_gnorm.
    p = dst / "gnorm.c"
    g = p.read_text(encoding="utf-8")
    if "#include <stdio.h>" not in g:
        g = g.replace("#include <stdint.h>",
                      "#include <stdint.h>\n#include <stdio.h>", 1)
    old_g = "  if (num128_lt(&tg, &norm)) {\n    return 1;\n  }"
    new_g = ("  fprintf(stderr, \"[gnorm] tg=%llu:%llu:%llu:%llu "
             "norm=%llu:%llu:%llu:%llu\\n\", "
             "(unsigned long long)tg.coeff[3], (unsigned long long)tg.coeff[2], "
             "(unsigned long long)tg.coeff[1], (unsigned long long)tg.coeff[0], "
             "(unsigned long long)norm.coeff[3], (unsigned long long)norm.coeff[2], "
             "(unsigned long long)norm.coeff[1], (unsigned long long)norm.coeff[0]);\n"
             + old_g)
    if old_g not in g:
        raise SystemExit("no encontre num128_lt en gnorm.c")
    p.write_text(g.replace(old_g, new_g, 1), encoding="utf-8")

    # 1) D = 128
    p = dst / "poly_param.h"
    t = p.read_text(encoding="utf-8")
    t, cnt = re.subn(r"#define\s+D\s+64LL", "#define D %dLL" % D, t, count=1)
    if cnt != 1:
        raise SystemExit("no pude fijar D en poly_param.h")
    p.write_text(t, encoding="utf-8")

    # M del esquema (rango funcional)
    p = dst / "param.h"
    t = p.read_text(encoding="utf-8")
    t, cnt = re.subn(r"(#define\s+M\s+)38LL", r"\g<1>%dLL" % M, t, count=1)
    if cnt != 1:
        raise SystemExit("no pude fijar M en param.h")
    p.write_text(t, encoding="utf-8")

    # 2) tamanos dependientes de D
    p = dst / "param.h"
    t = p.read_text(encoding="utf-8")
    sizes = {
        "REJ_UNIF_Q_G0_LEN": slen(M * N * D),
        "REJ_UNIF_Q_GR_LEN": slen((M + R) * N * D),
        "REJ_UNIF_Q_H_LEN": slen(M * D),
        "REJ_UNIF_QBIG_LEN": slen((M_HAT + V_HAT) * N_HAT * D),
        "REJ_UNIF_QBIG1_LEN": slen((N_HAT - 1) * D),
        "REJ_UNIF_QBIG2_LEN": slen((N_HAT - 1) * D),
        "SAMPLE_B_MD_LEN": slen(M * D),
        "SAMPLE_B_MHATD_LEN": slen(M_HAT * D),
        "SAMPLE_B_MAX_LEN": max(slen(M * D), slen(M_HAT * D)),
        "SAMPLE_BE_MHATD_LEN": slen(M_HAT * D),
        "SAMPLE_BE_V1D_LEN": slen(V_1 * D),
        "SAMPLE_BE_V2D_LEN": slen(V_2 * D),
        "SAMPLE_BE_MAX_LEN": max(slen(M_HAT * D), slen(V_1 * D), slen(V_2 * D)),
        "SAMPLE_BBIG_MD_LEN": slen(M * D),
        "SAMPLE_BBIG_MAX_LEN": max(slen(M * D), slen(M_HAT * D)),
        "SAMPLE_BBIG_K_PRIME_MD_LEN": slen(M * D),
        "SAMPLE_BBIG_HAT_MHATD_LEN": slen(M_HAT * D),
        "SAMPLE_BSTAR_MAX_LEN": slen(D),
        "SAMPLE_BA_LEN": slen(D),
        "SAMPLE_BR_LEN": slen(D),
        "SAMPLE_BAP_LEN": slen(D),
        "HASH_OUTPUT_LEN": 1024,
    }
    for name, val in sizes.items():
        t = setdef(t, name, val)

    # 3b) BR / BRP (opcional): BRP = BR - P se recalcula solo.
    if args.br_factor != 1:
        br = getdef(t, "BR") * args.br_factor
        t = setdef(t, "BR", br)
        # sample_bstar_d usa este factor con BARRETT_SHIFT_BSTAR = 24.
        t = setdef(t, "BARRETT_FACTOR_BR", (1 << 24) // (2 * br + 1))

    # 3) TG (opcional)
    if args.tg_factor != 1:
        tg = 0
        for i, name in enumerate(TG_NAMES):
            tg += getdef(t, name) << (32 * i)
        tg *= args.tg_factor
        for i, name in enumerate(TG_NAMES):
            t = setdef(t, name, (tg >> (32 * i)) & 0xFFFFFFFF)
    p.write_text(t, encoding="utf-8")

    # 4) mascara de parse_hash
    p = dst / "samplex.c"
    t = p.read_text(encoding="utf-8")
    if t.count("& 63") != 2:
        raise SystemExit("esperaba 2 usos de '& 63' en samplex.c")
    t = t.replace("& 63", "& (D - 1)")
    p.write_text(t, encoding="utf-8")

    # 5) multiplicacion de subanillo generica
    p = dst / "poly_q.c"
    t = p.read_text(encoding="utf-8")
    t, cnt = MULT_OLD_RE.subn(MULT_NEW, t, count=1)
    if cnt != 1:
        raise SystemExit("no encontre el bloque mult() original en poly_q.c")
    p.write_text(t, encoding="utf-8")

    # 5b) hash canónico: el challenge no debe depender de la representación
    # NTT (densa vs subanillos). Se hashean coeficientes.
    p = dst / "hash.c"
    t = p.read_text(encoding="utf-8")
    old = """  /* we convert each coordinates to a 32-bit integer within [0,q-1], where q can
   * be Q, QBIG1, or QBIG2 */
  for (i = 0; i < N_HAT; i++) {
    for (j = 0; j < D; j++) {
      store_32(hash_input_head, a[i].poly1[j]);
      store_32(hash_input_head + HASH_INPUT_BYTES, a[i].poly2[j]);"""
    new = """  /* Canonicalizacion: hashear coeficientes, no el dominio NTT, para que el
   * challenge sea independiente de la representacion (densa vs subanillos). */
  static POLY_QBIG a_can[N_HAT], b_can[N_HAT];
  for (i = 0; i < N_HAT; i++) {
    a_can[i] = a[i];
    b_can[i] = b[i];
    intt_qbig(a_can + i);
    intt_qbig(b_can + i);
  }
  for (i = 0; i < N_HAT; i++) {
    for (j = 0; j < D; j++) {
      store_32(hash_input_head, a_can[i].poly1[j]);
      store_32(hash_input_head + HASH_INPUT_BYTES, a_can[i].poly2[j]);"""
    if old not in t:
        raise SystemExit("no encontre el bloque de hash de a[]")
    t = t.replace(old, new, 1)
    t = t.replace("store_32(hash_input_head, b[i].poly1[j]);",
                  "store_32(hash_input_head, b_can[i].poly1[j]);", 1)
    t = t.replace("store_32(hash_input_head + HASH_INPUT_BYTES, b[i].poly2[j]);",
                  "store_32(hash_input_head + HASH_INPUT_BYTES, b_can[i].poly2[j]);", 1)
    p.write_text(t, encoding="utf-8")

    # 6) NTT de qbig: denso + subanillos (estilo Plus), con despacho por modo
    qbig_dense = (POLY_QBIG_C
                  .replace("void ntt_qbig(", "void ntt_qbig_dense(")
                  .replace("void intt_qbig(", "void intt_qbig_dense(")
                  .replace("void mult_rqbig(", "void mult_rqbig_dense("))
    (dst / "poly_qbig.c").write_text(qbig_dense, encoding="utf-8")
    (dst / "poly_qbig_sub.c").write_text(POLY_QBIG_SUB_C, encoding="utf-8")

    # 7) tablas del gadget (si se recuperaron)
    if args.gadget:
        entries = load_gadget(args.gadget)
        tg0 = []
        tg1 = []
        if D > 128:
            # modo bloques: raiz de orden 2*(D/4), existe en ambos modulos
            psi1 = root_2d(QBIG1, 2 * (D >> 2))
            psi2 = root_2d(QBIG2, 2 * (D >> 2))
            conv = ntt_py_sub
        else:
            psi1 = root_2d(QBIG1, D)
            psi2 = root_2d(QBIG2, D)
            conv = ntt_py
        for idx, (p1, p2) in enumerate(entries):
            e1 = conv(list(p1) + [0] * (D - D_OLD), QBIG1, psi1)
            e2 = conv(list(p2) + [0] * (D - D_OLD), QBIG2, psi2)
            if idx < N_SPENT_TABLES:
                tg0.append((e1, e2))
            else:
                tg1.append((e1, e2))
        p = dst / "gentdrow.c"
        t = p.read_text(encoding="utf-8")
        tg0 = tg0[:args.n_spent - 1]
        t = replace_table(
            t, "static const POLY_QBIG tg_0[",
            table_block("static const POLY_QBIG tg_0[N_SPENT - 1] = {", tg0))
        t = replace_table(
            t, "static const POLY_QBIG tg_1[",
            table_block("static const POLY_QBIG tg_1[S * R] = {", tg1))
        p.write_text(t, encoding="utf-8")
        print(f"OK: gadget regenerado ({len(tg0)} tg_0 + {len(tg1)} tg_1)")
    else:
        print("AVISO: sin --gadget; gentdrow.c queda con tablas D=64 "
              "(solo sirve para auditoria estatica, no para correr)")

    print(f"OK: {dst} generado para D={D} (tg_factor={args.tg_factor})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
