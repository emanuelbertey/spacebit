#!/usr/bin/env python3
"""Auditoria de seguridad de Au256 con sisrank=4, lwerank=5.
Corre los 3 algoritmos de ataque (variantes 0,1,2) como ASIS_est.sage y
valida los componentes de cada uno. Seguridad = min(coste) entre validos.
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
                       delta_BKZ, nvec_sieve, svp_quantum)

Q, D, S, LW, N = 2147221513, 256, 4, 5, 2
B = [63307776 - 448, 5229773 - 448, 98304 - 8, 10240 - 8, 10240]
M = [1, N + 1, S + LW, 1, N + 1]
MM = [D * x for x in M]
H = D * S


def main():
    sis.MIN_b, sis.STEPS_b, sis.STEPS_m = 40, 10, 10
    out = {"esquema": "Au256", "sisrank": S, "lwerank": LW, "n": N,
           "d": D, "h": H, "q": Q,
           "cotas_log2": [round(math.log2(x), 2) for x in B],
           "algoritmos": {}, "nota": "modelo SIS, 3 variantes de ataque"}

    for v in (0, 1, 2):
        p = sis.MSISParameterSet(D, sum(M), S, *B, *M, Q, norm="linf")
        w, b, c = SIS_optimize_attack(
            Q, D * sum(M), H, B[0], MM[0], B[1], MM[1], B[2], MM[2],
            B[3], MM[3], B[4], MM[4], cost_attack=SIS_linf_cost,
            cost_svp=svp_quantum, verbose=False, attack_variant=v)
        nq, n1 = w - H, H
        if v == 0:
            cc = [1.0] * 4
            i, j, L = construct_BKZ_shape(Q, nq, n1, b)
            glv = nq * math.log(Q)
        else:
            cc = [B[0] / B[t + 1] for t in range(4)] if v == 2 else [1.0] * 4
            i, j, L = construct_BKZ_shape_randomized(Q, nq, n1, b, *MM, *cc)
            dq = nq + n1
            d1 = min(MM[0], dq)
            d2 = min(max(0, dq - d1), MM[1])
            d3 = min(max(0, dq - d1 - d2), MM[2])
            d4 = min(max(0, dq - d1 - d2 - d3), MM[3])
            d5 = min(max(0, dq - d1 - d2 - d3 - d4), MM[4])
            glv = (nq * math.log(Q) + d2 * math.log(cc[0]) +
                   d3 * math.log(cc[1]) + d4 * math.log(cc[2]) +
                   d5 * math.log(cc[3]))
        l = math.exp(L[i]) if L else 0.0
        valido = l < Q
        deff = j - i + 1
        sigma = l / math.sqrt(deff) if deff else 0.0
        pendiente = -2 * math.log(delta_BKZ(b))
        tramo = [L[k + 1] - L[k] for k in range(i, min(j, len(L) - 1))]
        desv = max((abs(x - pendiente) for x in tramo), default=0.0)
        out["algoritmos"][f"variant_{v}"] = {
            "bkz": int(b), "dim_sis": int(w), "nq": nq, "dim_efectiva": deff,
            "det_ratio": round(sum(L) / glv, 8),
            "gso_desv": round(desv, 6),
            "rhf": round(float(delta_BKZ(b)), 6),
            "l": round(l, 2), "q": Q, "l_menor_q": int(valido),
            "c12_15": [round(x, 4) for x in cc],
            "cost_bits": round(float(c), 2), "aplicable": int(valido),
        }
    validos = [a["cost_bits"] for a in out["algoritmos"].values()
               if a["aplicable"]]
    out["coste_seguridad"] = min(validos) if validos else None
    out["algoritmo_ganador"] = [k for k, a in out["algoritmos"].items()
                                if a["aplicable"] and
                                a["cost_bits"] == out["coste_seguridad"]]

    dest = ROOT / "security" / "auditoria_seguridad_au256.json"
    dest.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
