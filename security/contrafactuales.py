#!/usr/bin/env python3
"""Contrafactuales en el MODELO (no tocan el protocolo).

Separan la causa de los 37.14 bits de Au256:
  - orden B4/B5 corregido
  - cotas simetricas
  - rango del modulo Au vs paper
  - cotas Au vs paper
Todos con busqueda gruesa (40/10/10) para que el optimizador pueda ver BKZ chico.
"""
import json
import math
import pathlib
import signal
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "security"))

import MSIS_security as sis
from MSIS_security import SIS_linf_cost, SIS_optimize_attack
from model_BKZ import (construct_BKZ_shape, construct_BKZ_shape_randomized,
                       delta_BKZ, svp_quantum)

D = 256
Q_AU, Q_PLUS = 2147221513, 3489661057
B_AU = [63307776 - 448, 5229773 - 448, 98304 - 8, 10240 - 8, 10240]
B_PLUS = [math.ceil(2 ** 30.3), math.ceil(2 ** 27.56), math.ceil(2 ** 22.64),
          math.ceil(2 ** 19.92), math.ceil(2 ** 16.6)]
M_AU, M_PLUS = [1, 3, 3, 1, 3], [1, 12, 10, 1, 12]


def esc(q, b, m, s, l):
    return {"q": q, "B": list(b), "m": list(m), "s": s, "l": l}


B4B5 = list(B_AU)
B4B5[3], B4B5[4] = B4B5[4], B4B5[3]
SIM = [B_AU[0]] * 5

ESCENARIOS = {
    "baseline_Au256": esc(Q_AU, B_AU, M_AU, 2, 1),
    "B4B5_ordenado": esc(Q_AU, B4B5, M_AU, 2, 1),
    "cotas_simetricas": esc(Q_AU, SIM, M_AU, 2, 1),
    "Au_cotas_rank_paper": esc(Q_AU, B_AU, M_PLUS, 5, 5),
    "paper_cotas_rank_Au": esc(Q_PLUS, B_PLUS, M_AU, 2, 1),
    "paper_completo": esc(Q_PLUS, B_PLUS, M_PLUS, 5, 5),
}


def geometria(p, v, w, b):
    q, s = p["q"], p["s"]
    B, m = p["B"], [D * x for x in p["m"]]
    h, nq, n1 = D * s, w - D * s, D * s
    if v == 0:
        c = [1.0] * 4
        i, j, L = construct_BKZ_shape(q, nq, n1, b)
        glv = nq * math.log(q)
    else:
        c = [B[0] / B[t + 1] for t in range(4)] if v == 2 else [1.0] * 4
        i, j, L = construct_BKZ_shape_randomized(q, nq, n1, b, *m, *c)
        dq = nq + n1
        d1 = min(m[0], dq)
        d2 = min(max(0, dq - d1), m[1])
        d3 = min(max(0, dq - d1 - d2), m[2])
        d4 = min(max(0, dq - d1 - d2 - d3), m[3])
        d5 = min(max(0, dq - d1 - d2 - d3 - d4), m[4])
        glv = (nq * math.log(q) + d2 * math.log(c[0]) + d3 * math.log(c[1]) +
               d4 * math.log(c[2]) + d5 * math.log(c[3]))
    l = math.exp(L[i]) if L else 0.0
    return {"bkz": int(b), "dim_sis": int(w), "h": h, "nq": nq,
            "glv": round(glv, 2), "ratio_det": round(sum(L) / glv, 6) if glv else None,
            "rhf": round(float(delta_BKZ(b)), 6), "l": round(l, 2),
            "l_menor_q": int(l < q), "c12_15": [round(x, 3) for x in c]}


class Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise Timeout()


def main():
    sis.MIN_b, sis.STEPS_b, sis.STEPS_m = 40, 10, 10
    signal.signal(signal.SIGALRM, _alarm)
    out = {"config": "gruesa 40/10/10", "nota": "min coste entre variantes con l<q",
           "escenarios": {}}
    for nombre, p in ESCENARIOS.items():
        q, B, m, s = p["q"], p["B"], p["m"], p["s"]
        mm = [D * x for x in m]
        h = D * s
        filas = []
        for v in (0, 1, 2):
            signal.alarm(90)
            try:
                params = sis.MSISParameterSet(D, sum(m), s, B[0], B[1], B[2],
                                              B[3], B[4], m[0], m[1], m[2],
                                              m[3], m[4], q, norm="linf")
                w, b, c = SIS_optimize_attack(
                    q, D * sum(m), h, B[0], mm[0], B[1], mm[1], B[2], mm[2],
                    B[3], mm[3], B[4], mm[4],
                    cost_attack=SIS_linf_cost, cost_svp=svp_quantum,
                    verbose=False, attack_variant=v)
                g = geometria(p, v, w, b)
                g.update(variante=v, cost_bits=round(float(c), 2),
                         valido=g["l_menor_q"])
            except Exception as exc:
                g = {"variante": v, "error": f"{type(exc).__name__}: {exc}"}
            finally:
                signal.alarm(0)
            filas.append(g)
        validos = [f["cost_bits"] for f in filas
                   if f.get("valido") and "cost_bits" in f]
        out["escenarios"][nombre] = {
            "log2_q": round(math.log2(q), 2),
            "cotas_log2": [round(math.log2(x), 2) for x in B],
            "m": m, "sisrank": s, "lwerank": p["l"],
            "corridas": filas,
            "coste_min_valido": min(validos) if validos else None,
        }
    dest = ROOT / "security" / "contrafactuales.json"
    dest.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
