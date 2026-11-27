#!/usr/bin/env python3
"""Separar cuanto viene de D y cuanto del bloque BKZ.

Fija BKZ y varia solo D (64 vs 256) con mismos ranks, cotas y m.
Luego barre ranks intermedios entre r3_l2 y r4_l3 buscando ~128 bits.
"""
import json
import math
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "security"))

import MSIS_security as sis
from MSIS_security import SIS_linf_cost, SIS_optimize_attack
from model_BKZ import construct_BKZ_shape, construct_BKZ_shape_randomized, delta_BKZ, svp_quantum

Q = 2147221513
B = [63307776 - 448, 5229773 - 448, 98304 - 8, 10240 - 8, 10240]


def correr(d, s, l, b_fijo=None):
    m = [1, 3, s + l, 1, 3]
    mm = [d * x for x in m]
    h = d * s
    filas = []
    for v in (0, 1, 2):
        p = sis.MSISParameterSet(d, sum(m), s, *B, *m, Q, norm="linf")
        if b_fijo is None:
            sis.MIN_b, sis.STEPS_b = 40, 10
        else:
            # Fijar b PERO dejar que el optimizador elija w.
            sis.MIN_b, sis.STEPS_b = b_fijo, 1
        w, b, c = SIS_optimize_attack(
            Q, d * sum(m), h, B[0], mm[0], B[1], mm[1], B[2], mm[2],
            B[3], mm[3], B[4], mm[4], cost_attack=SIS_linf_cost,
            cost_svp=svp_quantum, verbose=False, attack_variant=v)
        if v == 0:
            i, j, Lg = construct_BKZ_shape(Q, w - h, h, b)
        else:
            cc = [B[0] / B[t + 1] for t in range(4)] if v == 2 else [1.0] * 4
            i, j, Lg = construct_BKZ_shape_randomized(Q, w - h, h, b, *mm, *cc)
        lvec = math.exp(Lg[i]) if Lg else 0.0
        filas.append({"v": v, "bkz": int(b), "h": h, "l": round(lvec, 2),
                      "valido": int(lvec < Q), "cost": round(float(c), 2)})
    validos = [f["cost"] for f in filas if f["valido"]]
    return {"d": d, "sisrank": s, "lwerank": l, "h": h, "m": m,
            "bkz_fijo": b_fijo, "corridas": filas,
            "coste_min_valido": min(validos) if validos else None}


def main():
    sis.MIN_b, sis.STEPS_b, sis.STEPS_m = 40, 10, 10
    out = {"parte1_BKZ_fijo": {}, "parte2_ranks_intermedios": {}}
    # Parte 1: mismo BKZ, solo cambia D
    for b in (40, 100, 140, 160):
        out["parte1_BKZ_fijo"][f"bkz_{b}"] = {
            "Au64": correr(64, 2, 1, b_fijo=b),
            "Au256": correr(256, 2, 1, b_fijo=b),
        }
    # Parte 2: ranks intermedios buscando ~128 (sin optimizar, b=350 fijo)
    for s, l in ((3, 2), (3, 3), (4, 2), (4, 3)):
        out["parte2_ranks_intermedios"][f"r{s}l{l}"] = correr(256, s, l, b_fijo=350)
    p = ROOT / "security" / "aislar_d_vs_bkz.json"
    p.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
