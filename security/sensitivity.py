#!/usr/bin/env python3
"""Sensibilidad del estimador MSIS ante las cotas del Au256-baseline.

No cambia parametros del protocolo: solo varia una cota a la vez en el
modelo para identificar cual limita el coste del ataque (cost_pq).

Uso: python security/sensitivity.py
Salida: security/sensitivity_au256.json
"""
import json
import math
import pathlib
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import MSIS_security as sis  # noqa: E402
from MSIS_security import MSISParameterSet, MSIS_summarize_attacks  # noqa: E402

sis.MIN_b = 40
sis.STEPS_b = 10
sis.STEPS_m = 10

N, LWERANK, SISRANK = 2, 1, 2
Q, D = 2147221513, 256
B1, B2 = 63307776 - 448, 5229773 - 448
B3, B4, B5 = 98304 - 8, 10240 - 8, 10240
M1, M2, M3, M4, M5 = 1, 3, SISRANK + LWERANK, 1, 3
NOMBRES = ("ZB_NORM", "ZC_NORM", "BRP", "BAP", "BA")


def make_params(bounds):
    return MSISParameterSet(D, M1 + M2 + M3 + M4 + M5, SISRANK,
                            *bounds, M1, M2, M3, M4, M5, Q, norm="linf")


def cost(bounds, variant):
    m, b, c = MSIS_summarize_attacks(make_params(bounds), attack_variant=variant)
    return round(float(c), 2), int(b), int(m)


def main():
    base = (B1, B2, B3, B4, B5)
    t0 = time.time()
    report = {
        "baseline": "Au256-baseline",
        "nota": ("Solo el MODELO: se divide una cota por 2 a la vez. "
                 "No se cambian los parametros funcionales."),
        "d": D,
        "cotas_log2": [round(math.log2(x), 2) for x in base],
        "orden_esperado": "B1 >= B2 >= B3 >= B4 >= B5",
        "orden_cumplido": all(base[i] >= base[i + 1] for i in range(4)),
    }

    variantes = {}
    for variant in (1, 2):
        c0, b0, m0 = cost(base, variant)
        filas = []
        for i, nombre in enumerate(NOMBRES):
            nueva = list(base)
            nueva[i] = max(1, nueva[i] // 2)
            c1, b1, m1 = cost(tuple(nueva), variant)
            filas.append({
                "cota": nombre,
                "log2_original": round(math.log2(base[i]), 2),
                "log2_mitad": round(math.log2(nueva[i]), 2),
                "cost_original": c0,
                "cost_mitad": c1,
                "ganancia_bits": round(c1 - c0, 2),
                "bloque_original": b0,
                "bloque_mitad": b1,
            })
        variantes[f"variant_{variant}"] = {
            "cost_baseline": c0,
            "bloque_baseline": b0,
            "dim_sis_baseline": m0,
            "sensibilidad": filas,
        }
    report["variantes"] = variantes
    report["segundos"] = round(time.time() - t0, 1)

    out = HERE / "sensitivity_au256.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"\n-> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
