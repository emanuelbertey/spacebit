#!/usr/bin/env python3
"""Comparacion controlada A/B/C/D: separa efecto de parametros y de optimizador.

A = MatRiCT+ originales, busqueda fina (300/5/5)
B = Au256 actuales,      busqueda fina (300/5/5)
C = Au256 actuales,      busqueda gruesa (40/10/10)
D = MatRiCT+ originales, busqueda gruesa (40/10/10)

Por cada variante registra geometria: dim efectiva, det (glv vs perfil),
RHF, longitud l, modulo, condicion l<q y coste. El nivel de seguridad es el
MINIMO entre los ataques validos (l<q).
"""
import json
import math
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "security"))

import MSIS_security as sis
from MSIS_security import SIS_linf_cost, SIS_optimize_attack
from model_BKZ import (construct_BKZ_shape, construct_BKZ_shape_randomized,
                       delta_BKZ, svp_quantum)

# --- conjuntos de parametros ---
PLUS = {
    "nombre": "MatRiCT+ (paper)",
    "q": 3489661057, "d": 256, "sisrank": 5, "lwerank": 5, "n": 11,
    "B": [math.ceil(2 ** 30.3), math.ceil(2 ** 27.56), math.ceil(2 ** 22.64),
          math.ceil(2 ** 19.92), math.ceil(2 ** 16.6)],
    "m": [1, 12, 10, 1, 12],
}
AU256 = {
    "nombre": "Au256-baseline",
    "q": 2147221513, "d": 256, "sisrank": 2, "lwerank": 1, "n": 2,
    "B": [63307776 - 448, 5229773 - 448, 98304 - 8, 10240 - 8, 10240],
    "m": [1, 3, 3, 1, 3],
}
CONFIGS = {"fina": (300, 5, 5), "gruesa": (40, 10, 10)}
EXPERIMENTOS = [("A", PLUS, "fina"), ("B", AU256, "fina"),
                ("C", AU256, "gruesa"), ("D", PLUS, "gruesa")]


def geometria(p, variante, w, b):
    q, d, s = p["q"], p["d"], p["sisrank"]
    B, m = p["B"], [d * x for x in p["m"]]
    h = d * s
    nq, n1 = w - h, h
    if variante == 0:
        c = [1.0] * 4
        i, j, L = construct_BKZ_shape(q, nq, n1, b)
        glv = nq * math.log(q)
    else:
        c = [B[0] / B[t + 1] for t in range(4)] if variante == 2 else [1.0] * 4
        i, j, L = construct_BKZ_shape_randomized(q, nq, n1, b, *m, *c)
        dq = nq + n1
        d1 = min(m[0], dq)
        d2 = min(max(0, dq - d1), m[1])
        d3 = min(max(0, dq - d1 - d2), m[2])
        d4 = min(max(0, dq - d1 - d2 - d3), m[3])
        d5 = min(max(0, dq - d1 - d2 - d3 - d4), m[4])
        glv = (nq * math.log(q) + d2 * math.log(c[0]) + d3 * math.log(c[1]) +
               d4 * math.log(c[2]) + d5 * math.log(c[3]))
    lv = sum(L)
    l = math.exp(L[i]) if L else 0.0
    return {"bkz": int(b), "dim_sis": int(w), "nq": nq, "n1": n1,
            "dim_efectiva": int(j - i + 1),
            "glv": round(glv, 3), "sum_L": round(lv, 3),
            "ratio_det": round(lv / glv, 6) if glv else None,
            "rhf": round(float(delta_BKZ(b)), 6),
            "l": round(l, 2), "q": q, "l_menor_q": int(l < q),
            "c12_15": [round(x, 4) for x in c]}


def main():
    out = {"nota": "Seguridad = min(coste) entre variantes con l<q",
           "experimentos": {}}
    for nombre, p, cfg in EXPERIMENTOS:
        minb, sb, sm = CONFIGS[cfg]
        sis.MIN_b, sis.STEPS_b, sis.STEPS_m = minb, sb, sm
        q, d, s = p["q"], p["d"], p["sisrank"]
        B, m = p["B"], p["m"]
        h = d * s
        mmax = d * sum(m)
        mm = [d * x for x in m]
        filas = []
        for v in (0, 1, 2):
            try:
                params = sis.MSISParameterSet(d, sum(m), s, B[0], B[1], B[2],
                                              B[3], B[4], m[0], m[1], m[2],
                                              m[3], m[4], q, norm="linf")
                w, b, c = SIS_optimize_attack(
                    q, mmax, h, B[0], mm[0], B[1], mm[1], B[2], mm[2],
                    B[3], mm[3], B[4], mm[4],
                    cost_attack=SIS_linf_cost, cost_svp=svp_quantum,
                    verbose=False, attack_variant=v)
                g = geometria(p, v, w, b)
                g.update(variante=v, cost_bits=round(float(c), 2), valido=int(g["l_menor_q"]))
            except Exception as exc:
                g = {"variante": v, "error": f"{type(exc).__name__}: {exc}"}
            filas.append(g)
        validos = [f["cost_bits"] for f in filas
                   if f.get("valido") and "cost_bits" in f]
        out["experimentos"][nombre] = {
            "parametros": p["nombre"], "config": cfg,
            "search": {"MIN_b": minb, "STEPS_b": sb, "STEPS_m": sm},
            "log2_q": round(math.log2(q), 2),
            "cotas_log2": [round(math.log2(x), 2) for x in B],
            "corridas": filas,
            "coste_min_valido": min(validos) if validos else None,
        }

    dest = ROOT / "security" / "compare_abcd.json"
    dest.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
