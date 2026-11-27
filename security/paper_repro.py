#!/usr/bin/env python3
"""Reproduce ASIS_est.sage del paper MatRiCT+ con el estimador ORIGINAL y
compararlo con nuestro estimador modificado (security/).

Parametros publicados (ASIS_est.sage):
  N=11, lwerank=5, sisrank=5, d=256,
  q = 3489661057  (es qhat, el paper Also usa q=167770241),
  B1..B5 = ceil(2^30.3), ceil(2^27.56), ceil(2^22.64),
           ceil(2^19.92), ceil(2^16.6),
  m1=1, m2=N+1, m3=sisrank+lwerank, m4=1, m5=N+1.

Uso: python security/paper_repro.py
Salida: security/paper_repro.json
"""
import importlib.util
import json
import math
import pathlib
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent.parent
ORIG = ROOT / "main" / "repos_originales" / "matrict_plus" / "asymmetric_sis_estimate"
OURS = ROOT / "security"


def cargar(nombre, directorio):
    sys.path.insert(0, str(directorio))
    spec = importlib.util.spec_from_file_location(
        nombre, directorio / "MSIS_security.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


N, LWERANK, SISRANK, D = 11, 5, 5, 256
B = [math.ceil(2 ** 30.3), math.ceil(2 ** 27.56), math.ceil(2 ** 22.64),
     math.ceil(2 ** 19.92), math.ceil(2 ** 16.6)]
M = [1, N + 1, SISRANK + LWERANK, 1, N + 1]
MODULOS = {"qhat_3489661057": 3489661057, "q_167770241": 167770241}


def correr(mod, q, variante, config):
    mod.MIN_b, mod.STEPS_b, mod.STEPS_m = config
    params = mod.MSISParameterSet(D, sum(M), SISRANK, B[0], B[1], B[2],
                                  B[3], B[4], M[0], M[1], M[2], M[3], M[4],
                                  q, norm="linf")
    m, b, c = mod.MSIS_summarize_attacks(params, attack_variant=variante)
    from model_BKZ import delta_BKZ
    return {"dim_sis": int(m), "bkz": int(b),
            "rhf": round(float(delta_BKZ(b)), 5),
            "cost_bits": round(float(c), 2)}


def main():
    t0 = time.time()
    report = {
        "fuente": "ASIS_est.sage (MatRiCT+, IEEE S&P 2022)",
        "parametros": {"N": N, "lwerank": LWERANK, "sisrank": SISRANK,
                       "d": D, "m": M,
                       "B_log2": [round(math.log2(x), 2) for x in B]},
        "rhf_paper": "1.0043 - 1.00469 (M-SIS/M-LWE, segun el paper)",
        "objetivo_paper": "~128 bits poscuanticos",
        "corridas": [],
    }
    configs = {"paper": (300, 5, 5), "nuestro_actual": (40, 10, 10)}
    report["configs"] = {k: {"MIN_b": v[0], "STEPS_b": v[1], "STEPS_m": v[2]}
                         for k, v in configs.items()}

    orig = cargar("orig_MSIS", ORIG)
    ours = cargar("our_MSIS", OURS)
    for nombre_mod, q in MODULOS.items():
        for nombre_est, mod in (("original_paper", orig), ("nuestro_modificado", ours)):
            for nombre_cfg, cfg in configs.items():
                for variante in (0, 1, 2):
                    r = correr(mod, q, variante, cfg)
                    r.update(modulo=nombre_mod, estimador=nombre_est,
                             config=nombre_cfg, variante=variante)
                    report["corridas"].append(r)

    report["segundos"] = round(time.time() - t0, 1)
    out = ROOT / "security" / "paper_repro.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"\n-> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
