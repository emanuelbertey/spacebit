#!/usr/bin/env python3
"""Comparacion correcta: MatRiCT-Au D=64 vs el port Au256.

Mismo esquema (sisrank=2, lwerank=1, n=2), mismas cotas, mismo m,
mismo optimizador y mismo criterio de validez (l<q). Solo cambia D.
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

Q = 2147221513
B = [63307776 - 448, 5229773 - 448, 98304 - 8, 10240 - 8, 10240]
M = [1, 3, 3, 1, 3]
SISRANK = 2

CASOS = {"Au64_original": 64, "Au256_port": 256}


def geometria(d, v, w, b):
    m = [d * x for x in M]
    h = d * SISRANK
    nq, n1 = w - h, h
    if v == 0:
        c = [1.0] * 4
        i, j, L = construct_BKZ_shape(Q, nq, n1, b)
        glv = nq * math.log(Q)
    else:
        c = [B[0] / B[t + 1] for t in range(4)] if v == 2 else [1.0] * 4
        i, j, L = construct_BKZ_shape_randomized(Q, nq, n1, b, *m, *c)
        dq = nq + n1
        d1 = min(m[0], dq)
        d2 = min(max(0, dq - d1), m[1])
        d3 = min(max(0, dq - d1 - d2), m[2])
        d4 = min(max(0, dq - d1 - d2 - d3), m[3])
        d5 = min(max(0, dq - d1 - d2 - d3 - d4), m[4])
        glv = (nq * math.log(Q) + d2 * math.log(c[0]) + d3 * math.log(c[1]) +
               d4 * math.log(c[2]) + d5 * math.log(c[3]))
    l = math.exp(L[i]) if L else 0.0
    return {"bkz": int(b), "dim_sis": int(w), "d": d, "h": h, "nq": nq,
            "dim_efectiva": int(j - i + 1),
            "ratio_det": round(sum(L) / glv, 6) if glv else None,
            "rhf": round(float(delta_BKZ(b)), 6),
            "l": round(l, 2), "q": Q, "l_menor_q": int(l < Q),
            "c12_15": [round(x, 3) for x in c]}


def main():
    sis.MIN_b, sis.STEPS_b, sis.STEPS_m = 40, 10, 10
    out = {"nota": "mismo esquema, mismas cotas, solo cambia D",
           "config": "gruesa 40/10/10",
           "cotas_log2": [round(math.log2(x), 2) for x in B],
           "casos": {}}
    for nombre, d in CASOS.items():
        m = [d * x for x in M]
        h = d * SISRANK
        filas = []
        for v in (0, 1, 2):
            params = sis.MSISParameterSet(d, sum(M), SISRANK, B[0], B[1], B[2],
                                          B[3], B[4], M[0], M[1], M[2], M[3],
                                          M[4], Q, norm="linf")
            w, b, c = SIS_optimize_attack(
                Q, d * sum(M), h, B[0], m[0], B[1], m[1], B[2], m[2],
                B[3], m[3], B[4], m[4],
                cost_attack=SIS_linf_cost, cost_svp=svp_quantum,
                verbose=False, attack_variant=v)
            g = geometria(d, v, w, b)
            g.update(variante=v, cost_bits=round(float(c), 2), valido=g["l_menor_q"])
            filas.append(g)
        validos = [f["cost_bits"] for f in filas if f["valido"]]
        out["casos"][nombre] = {"d": d, "h": h, "corridas": filas,
                                "coste_min_valido": min(validos) if validos else None}

    dest = ROOT / "security" / "au64_vs_au256.json"
    dest.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
