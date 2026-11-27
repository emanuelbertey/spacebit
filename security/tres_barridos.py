#!/usr/bin/env python3
"""Tres barridos: ranks por separado, n, y cotas. Validez: l<q."""
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

Q, D = 2147221513, 256
B0 = [63307776 - 448, 5229773 - 448, 98304 - 8, 10240 - 8, 10240]


def caso(nombre, s, lwer, n, B):
    m = [1, n + 1, s + lwer, 1, n + 1]
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
        lv = math.exp(Lg[i]) if Lg else 0.0
        filas.append({"v": v, "bkz": int(b), "h": h,
                      "l": round(lv, 2), "valido_l_menor_q": int(lv < Q),
                      "cost": round(float(c), 2)})
    validos = [f["cost"] for f in filas if f["valido_l_menor_q"]]
    return {"nombre": nombre, "sisrank": s, "lwerank": lwer, "n": n,
            "h": h, "m": m, "cotas_log2": [round(math.log2(x), 2) for x in B],
            "corridas": filas, "coste_min_valido": min(validos) if validos else None}


def main():
    sis.MIN_b, sis.STEPS_b, sis.STEPS_m = 40, 10, 10
    out = {"d": D, "validez": "l<q", "experimentos": {}}

    out["experimentos"]["1_sisrank_fijo_lwerank_varia"] = [
        caso(f"s2_l{l}_n2", 2, l, 2, B0) for l in (1, 2, 3)]
    out["experimentos"]["2_lwerank_fijo_sisrank_varia"] = [
        caso(f"s{s}_l1_n2", s, 1, 2, B0) for s in (2, 3, 4)]
    out["experimentos"]["3_n_varia"] = [
        caso(f"s2_l1_n{n}", 2, 1, n, B0) for n in (2, 3, 4, 5, 6, 8, 11)]
    out["experimentos"]["4_cotas_varia"] = [
        caso(f"cotas_x{f}_s2l1n2", 2, 1, 2,
             [max(1, int(x * f)) for x in B0]) for f in (0.5, 0.75, 1.0, 1.25)]

    p = ROOT / "security" / "tres_barridos.json"
    p.write_text(json.dumps(out, indent=2), encoding="utf-8")
    resumen = {k: [(x["nombre"], x["coste_min_valido"]) for x in v]
               for k, v in out["experimentos"].items()}
    print(json.dumps(resumen, indent=2))
    print(f"-> {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
