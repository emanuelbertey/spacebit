#!/usr/bin/env python3
"""OFICIAL (MatRiCT-Au D=64) vs NUESTRO (port D=256).
Mismo esquema, mismas cotas, mismo optimizador, misma validez (l<q).
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
S, L = 2, 1

CASOS = [("OFICIAL_Au_D64", 64), ("NUESTRO_Au_D256", 256)]


def main():
    sis.MIN_b, sis.STEPS_b, sis.STEPS_m = 40, 10, 10
    out = {"comparacion": "OFICIAL MatRiCT-Au D=64 vs NUESTRO port D=256",
           "mismo_esquema": {"sisrank": S, "lwerank": L, "n": 2},
           "mismas_cotas_log2": [round(math.log2(x), 2) for x in B],
           "config": "gruesa 40/10/10", "casos": {}}
    for nombre, d in CASOS:
        mm = [d * x for x in M]
        h = d * S
        filas = []
        for v in (0, 1, 2):
            p = sis.MSISParameterSet(d, sum(M), S, *B, *M, Q, norm="linf")
            w, b, c = SIS_optimize_attack(
                Q, d * sum(M), h, B[0], mm[0], B[1], mm[1], B[2], mm[2],
                B[3], mm[3], B[4], mm[4], cost_attack=SIS_linf_cost,
                cost_svp=svp_quantum, verbose=False, attack_variant=v)
            if v == 0:
                i, j, Lg = construct_BKZ_shape(Q, w - h, h, b)
            else:
                cc = [B[0] / B[t + 1] for t in range(4)] if v == 2 else [1.0] * 4
                i, j, Lg = construct_BKZ_shape_randomized(Q, w - h, h, b, *mm, *cc)
            lvec = math.exp(Lg[i]) if Lg else 0.0
            filas.append({"variante": v, "bkz": int(b), "h": h,
                          "rhf": round(float(delta_BKZ(b)), 6),
                          "l": round(lvec, 2), "valido_l_menor_q": int(lvec < Q),
                          "cost_bits": round(float(c), 2)})
        validos = [f["cost_bits"] for f in filas if f["valido_l_menor_q"]]
        out["casos"][nombre] = {"d": d, "h": h, "corridas": filas,
                                "coste_min_valido": min(validos) if validos else None}
    a = out["casos"]["OFICIAL_Au_D64"]["coste_min_valido"]
    b = out["casos"]["NUESTRO_Au_D256"]["coste_min_valido"]
    out["resultado"] = ("NUESTRO queda ARRIBA" if b > a else "NUESTRO queda ABAJO")
    out["diferencia_bits"] = round(b - a, 2) if a is not None and b is not None else None
    p = ROOT / "security" / "test_oficial_vs_nuestro.json"
    p.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
