#!/usr/bin/env python3
"""Estimacion de seguridad MSIS con los scripts de main/ (Python puro).

Uso: python security/estimate.py [plus|au]
  plus: parametros publicados MatRiCT+ (ASIS_est.sage de main/)
  au:   parametros derivados de los topes del codigo Au (aproximado)

Emite security/estimate_<scheme>.json. Cada estimacion tiene un limite de
tiempo para no pasar del presupuesto del job (60 min).
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
from MSIS_security import MSISParameterSet, MSIS_summarize_attacks  # noqa: E402
from model_BKZ import delta_BKZ  # noqa: E402

# El paper fija MIN_b=300 para sus parametros; con d=64/256 hay menos
# bloques y el optimizador se queja. Se baja al minimo seguro.
sis.MIN_b = 60
sis.STEPS_b = 25
sis.STEPS_m = 25


class Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise Timeout()


def estimate(params, seconds=1200):
    tiene_alarma = hasattr(signal, "SIGALRM")
    if tiene_alarma:
        signal.signal(signal.SIGALRM, _alarm)
        signal.alarm(seconds)
    t0 = time.time()
    try:
        out = []
        for variant in (0, 1, 2):
            m, b, c = MSIS_summarize_attacks(params, attack_variant=variant)
            out.append({"variant": variant, "dim": m, "bkz_block": b,
                        "rhf": round(delta_BKZ(b), 5), "cost_bits": round(c, 2)})
        if tiene_alarma:
            signal.alarm(0)
        return out, time.time() - t0, None
    except Timeout:
        if tiene_alarma:
            signal.alarm(0)
        return [], time.time() - t0, f"timeout tras {seconds}s"


def main():
    scheme = sys.argv[1] if len(sys.argv) > 1 else "plus"
    if scheme == "plus":
        # Parametros publicados (ASIS_est.sage, main/repos_originales)
        N, lwerank, sisrank = 11, 5, 5
        q, d = 3489661057, 256
        B1, B2, B3, B4, B5 = (int(2**30.3), int(2**27.56), int(2**22.64),
                              int(2**19.92), int(2**16.6))
        m1, m2, m3, m4, m5 = 1, N + 1, sisrank + lwerank, 1, N + 1
        nota = "parametros publicados MatRiCT+ (paper IEEE S&P 2022)"
    elif scheme in ("au", "au256"):
        # Derivado de auditor/param.h: B1=ZB_NORM, B2=ZC_NORM, B3=BRP,
        # B4=BAP, B5=BA (cotas del verificador). Aproximado, no del paper.
        # au256 usa las MISMAS cotas pero con D=256 (anillo de auditor
        # portado a bloques). No es una estimacion oficial ni una prueba.
        N, lwerank, sisrank = 2, 1, 2
        q, d = 2147221513, (64 if scheme == "au" else 256)
        B1, B2 = 63307776 - 448, 5229773 - 448
        B3, B4, B5 = 98304 - 8, 10240 - 8, 10240
        m1, m2, m3, m4, m5 = 1, 3, sisrank + lwerank, 1, 3
        nota = ("APROXIMADO: cotas del codigo Au con D=256 (bloques)"
                if scheme == "au256" else
                "APROXIMADO: cotas del codigo Au, no estimacion oficial del paper")
    else:
        print("uso: estimate.py [plus|au]")
        return 2

    if scheme.startswith("au"):
        sis.MIN_b = 40
        sis.STEPS_b = 10
        sis.STEPS_m = 10

    params = MSISParameterSet(d, m1 + m2 + m3 + m4 + m5, sisrank,
                              B1, B2, B3, B4, B5, m1, m2, m3, m4, m5,
                              q, norm="linf")
    results, elapsed, error = estimate(params, seconds=900)
    report = {"scheme": scheme, "baseline": ("Au256-baseline" if scheme == "au256"
                                             else None),
              "nota": nota, "log2_q": round(math.log2(q), 2),
              "d": d, "bounds_log2": [round(math.log2(b), 2)
                                      for b in (B1, B2, B3, B4, B5)],
              "results": results, "segundos": round(elapsed, 1), "error": error}
    out = HERE / f"estimate_{scheme}.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"\n-> {out}")
    return 1 if error else 0


if __name__ == "__main__":
    sys.exit(main())
