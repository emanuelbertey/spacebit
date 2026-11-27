#!/usr/bin/env python3
"""Auditoria del estimador MSIS del Au256-baseline.

No cambia parametros del protocolo. Verifica la suposicion de monotonía de
cotas, corre los tres modelos de SVP (cuantico/clasico/plausible), las tres
variantes de ataque y compara la busqueda actual (gruesa) con la del paper.

Uso: python security/estimate_audit.py
Salida: security/estimate_audit.json
"""
import json
import math
import pathlib
import signal
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import MSIS_security as sis  # noqa: E402
from MSIS_security import SIS_linf_cost, SIS_optimize_attack  # noqa: E402
from model_BKZ import svp_classical, svp_quantum, svp_plausible  # noqa: E402

D, Q, SISRANK, LWERANK = 256, 2147221513, 2, 1
B = (63307776 - 448, 5229773 - 448, 98304 - 8, 10240 - 8, 10240)
M = (1, 3, SISRANK + LWERANK, 1, 3)
NOMBRES = ("ZB_NORM", "ZC_NORM", "BRP", "BAP", "BA")

H = D * SISRANK
MAX_W = D * sum(M)
MM = [D * m for m in M]

MODELOS = {"quantum": svp_quantum, "classic": svp_classical,
           "plausible": svp_plausible}
CONFIGS = {"actual": (40, 10, 10), "paper": (300, 5, 5)}


class Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise Timeout()


def correr(bounds, min_b, steps_b, steps_m, modelo, variante, segundos=600):
    sis.MIN_b, sis.STEPS_b, sis.STEPS_m = min_b, steps_b, steps_m
    tiene = hasattr(signal, "SIGALRM")
    if tiene:
        signal.signal(signal.SIGALRM, _alarm)
        signal.alarm(segundos)
    t0 = time.time()
    try:
        w, b, c = SIS_optimize_attack(
            Q, MAX_W, H,
            bounds[0], MM[0], bounds[1], MM[1], bounds[2], MM[2],
            bounds[3], MM[3], bounds[4], MM[4],
            cost_attack=SIS_linf_cost, cost_svp=MODELOS[modelo],
            verbose=False, attack_variant=variante)
        if tiene:
            signal.alarm(0)
        return {"dim_sis": int(w), "bkz": int(b), "cost_bits": round(float(c), 2),
                "segundos": round(time.time() - t0, 1), "error": None}
    except Timeout:
        if tiene:
            signal.alarm(0)
        return {"dim_sis": None, "bkz": None, "cost_bits": None,
                "segundos": round(time.time() - t0, 1),
                "error": f"timeout tras {segundos}s"}
    except Exception as exc:
        if tiene:
            signal.alarm(0)
        return {"dim_sis": None, "bkz": None, "cost_bits": None,
                "segundos": round(time.time() - t0, 1),
                "error": f"{type(exc).__name__}: {exc}"}


def main():
    t0 = time.time()
    monotonico = all(B[i] >= B[i + 1] for i in range(4))
    report = {
        "baseline": "Au256-baseline",
        "d": D,
        "q": Q,
        "h": H,
        "max_w": MAX_W,
        "cotas": [{"nombre": n, "valor": v, "log2": round(math.log2(v), 2),
                   "m_polinos": m}
                  for n, v, m in zip(NOMBRES, B, M)],
        "orden_B1_ge_B2_ge_B3_ge_B4_ge_B5": monotonico,
        "nota_orden": ("B4=BAP=10232 < B5=BA=10240: el modelo asume lo contrario. "
                       "No se reordenan valores artificialmente."),
        "configs": {k: {"MIN_b": v[0], "STEPS_b": v[1], "STEPS_m": v[2]}
                    for k, v in CONFIGS.items()},
        "corridas": [],
    }

    # Configuracion actual: todos los modelos y variantes.
    for variante in (0, 1, 2):
        for modelo in MODELOS:
            r = correr(B, *CONFIGS["actual"], modelo, variante)
            r.update(config="actual", modelo=modelo, variante=variante)
            report["corridas"].append(r)

    # Configuracion del paper: solo el ataque mas fuerte (variante 2).
    for modelo in ("quantum", "classic"):
        r = correr(B, *CONFIGS["paper"], modelo, 2, segundos=900)
        r.update(config="paper", modelo=modelo, variante=2)
        report["corridas"].append(r)

    report["segundos"] = round(time.time() - t0, 1)
    out = HERE / "estimate_audit.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"\n-> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
