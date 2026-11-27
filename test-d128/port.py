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

  /* inverse Cooley-Tukey */
  for (len = 2; len <= n; len <<= 1) {
    wl = modpow_q(winv, n / len, q);
    for (i = 0; i < n; i += len) {
      ww = 1;
      for (j = 0; j < (len >> 1); j++) {
        u = a[i + j];
        v = a[i + j + (len >> 1)];
        a[i + j] = (u + v) % q;
        a[i + j + (len >> 1)] = ((u + q - v) % q) * ww % q;
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
    lines = text.splitlines()
    start = None
    for i, line in enumerate(lines):
        if line.startswith(decl_start):
            start = i
            break
    if start is None:
        raise SystemExit(f"no encontre la tabla {decl_start}")
    end = None
    for i in range(start, len(lines)):
        if lines[i].strip() == "};":
            end = i
            break
    if end is None:
        raise SystemExit(f"no encontre el cierre de {decl_start}")
    return "\n".join(lines[:start] + new_block.splitlines() + lines[end + 1:]) + "\n"


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="auditor_d128")
    ap.add_argument("--gadget", default=None,
                    help="gadget_polys.txt recuperado con el codigo original D=64")
    ap.add_argument("--tg-factor", type=int, default=1)
    args = ap.parse_args()

    dst = ROOT / args.out
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(SRC, dst)

    # 1) D = 128
    p = dst / "poly_param.h"
    t = p.read_text(encoding="utf-8")
    t, cnt = re.subn(r"#define\s+D\s+64LL", "#define D %dLL" % D, t, count=1)
    if cnt != 1:
        raise SystemExit("no pude fijar D=128 en poly_param.h")
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

    # 6) NTT de qbig generico
    (dst / "poly_qbig.c").write_text(POLY_QBIG_C, encoding="utf-8")

    # 7) tablas del gadget (si se recuperaron)
    if args.gadget:
        entries = load_gadget(args.gadget)
        psi1 = root_2d(QBIG1, D)
        psi2 = root_2d(QBIG2, D)
        tg0 = []
        tg1 = []
        for idx, (p1, p2) in enumerate(entries):
            e1 = ntt_py(list(p1) + [0] * (D - D_OLD), QBIG1, psi1)
            e2 = ntt_py(list(p2) + [0] * (D - D_OLD), QBIG2, psi2)
            if idx < N_SPENT_TABLES:
                tg0.append((e1, e2))
            else:
                tg1.append((e1, e2))
        p = dst / "gentdrow.c"
        t = p.read_text(encoding="utf-8")
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
