#!/usr/bin/env python3
"""Reproduce ASIS_est.sage con el estimador ORIGINAL del paper.
Simple: un solo estimador, parametros publicados, 3 variantes, 2 moduli.
"""
import json
import math
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
ORIG = ROOT / "main" / "repos_originales" / "matrict_plus" / "asymmetric_sis_estimate"
sys.path.insert(0, str(ORIG))

import MSIS_security as sis
from MSIS_security import MSISParameterSet, MSIS_summarize_attacks
from model_BKZ import delta_BKZ

N, LWERANK, SISRANK, D = 11, 5, 5, 256
B = [math.ceil(2 ** 30.3), math.ceil(2 ** 27.56), math.ceil(2 ** 22.64),
     math.ceil(2 ** 19.92), math.ceil(2 ** 16.6)]
M = [1, N + 1, SISRANK + LWERANK, 1, N + 1]
MODULOS = {"qhat": 3489661057, "q": 167770241}

sis.MIN_b, sis.STEPS_b, sis.STEPS_m = 300, 5, 5

out = {"fuente": "ASIS_est.sage", "d": D, "modulos": {}}
for nombre, q in MODULOS.items():
    params = MSISParameterSet(D, sum(M), SISRANK, B[0], B[1], B[2], B[3], B[4],
                              M[0], M[1], M[2], M[3], M[4], q, norm="linf")
    filas = []
    for v in (0, 1, 2):
        try:
            m, b, c = MSIS_summarize_attacks(params, attack_variant=v)
            filas.append({"variante": v, "dim_sis": int(m), "bkz": int(b),
                          "rhf": round(float(delta_BKZ(b)), 5),
                          "cost_bits": round(float(c), 2)})
        except Exception as exc:
            filas.append({"variante": v,
                          "error": f"{type(exc).__name__}: {exc}"})
    out["modulos"][nombre] = {"q": q,
                              "log2_q": round(math.log2(q), 2),
                              "resultados": filas}

p = ROOT / "security" / "paper_repro.json"
p.write_text(json.dumps(out, indent=2), encoding="utf-8")
print(json.dumps(out, indent=2))
