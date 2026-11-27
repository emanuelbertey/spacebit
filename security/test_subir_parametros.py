#!/usr/bin/env python3
"""Subir parametros del modulo (ranks) con D=256 fijo: cuanto suben los bits."""
import json
import math
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "security"))

import MSIS_security as sis
from MSIS_security import SIS_linf_cost, SIS_optimize_attack
from model_BKZ import construct_BKZ_shape, construct_BKZ_shape_randomized, delta_BKZ, svp_quantum

Q, D = 2147221513, 256
B = [63307776 - 448, 5229773 - 448, 98304 - 8, 10240 - 8, 10240]
CASOS = [("r2_l1_actual", 2, 1), ("r3_l2", 3, 2), ("r4_l3", 4, 3), ("r5_l5", 5, 5)]


def main():
    sis.MIN_b, sis.STEPS_b, sis.STEPS_m = 40, 10, 10
    out = {"d": D, "mismas_cotas": [round(math.log2(x), 2) for x in B],
           "config": "gruesa 40/10/10", "casos": {}}
    for nombre, s, l in CASOS:
        m = [1, 3, s + l, 1, 3]
        mm = [D * x for x in m]
        h = D * s
        filas = []
        for v in (0, 1, 2):
            p = sis.MSISParameterSet(D, sum(m), s, *B, *m, Q, norm="linf")
            w, b, c = SIS_optimize_attack(
                Q, D * sum(m), h, B[0], mm[0], B[1], mm[1], B[2], mm[2],
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
                          "l": round(lvec, 2), "valido": int(lvec < Q),
                          "cost_bits": round(float(c), 2)})
        validos = [f["cost_bits"] for f in filas if f["valido"]]
        out["casos"][nombre] = {"sisrank": s, "lwerank": l, "h": h, "m": m,
                                "corridas": filas,
                                "coste_min_valido": min(validos) if validos else None}
    p = ROOT / "security" / "test_subir_parametros.json"
    p.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
