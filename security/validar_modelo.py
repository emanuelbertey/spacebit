#!/usr/bin/env python3
"""Bateria de validacion independiente de los componentes del estimador.

Instancia: MatRiCT+ publicado (ASIS_est.sage), optimos de paper_repro.json.
Reconstruye determinante, GSO, l, p_head, p_middle, epsilon, log2_R y coste
sin llamar a SIS_linf_cost, y compara contra el estimador.
"""
import json
import math
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "security"))

from model_BKZ import (construct_BKZ_shape, construct_BKZ_shape_randomized,
                       delta_BKZ, nvec_sieve, svp_quantum)
import MSIS_security as sis

Q, D, S = 3489661057, 256, 5
B = [math.ceil(2 ** 30.3), math.ceil(2 ** 27.56), math.ceil(2 ** 22.64),
     math.ceil(2 ** 19.92), math.ceil(2 ** 16.6)]
M = [1, 12, 10, 1, 12]
MM = [D * x for x in M]
H = D * S
OPTIMOS = {0: (2331, 315), 1: (2391, 350), 2: (2456, 340)}


def check(nombre, ok, detalle):
    return {"check": nombre, "ok": bool(ok), "detalle": detalle}


def validar(v, w, b):
    nq, n1 = w - H, H
    if v == 0:
        c = [1.0] * 4
        i, j, L = construct_BKZ_shape(Q, nq, n1, b)
        glv = nq * math.log(Q)
    else:
        c = [B[0] / B[t + 1] for t in range(4)] if v == 2 else [1.0] * 4
        i, j, L = construct_BKZ_shape_randomized(Q, nq, n1, b, *MM, *c)
        dq = nq + n1
        d1 = min(MM[0], dq)
        d2 = min(max(0, dq - d1), MM[1])
        d3 = min(max(0, dq - d1 - d2), MM[2])
        d4 = min(max(0, dq - d1 - d2 - d3), MM[3])
        d5 = min(max(0, dq - d1 - d2 - d3 - d4), MM[4])
        glv = (nq * math.log(Q) + d2 * math.log(c[0]) + d3 * math.log(c[1]) +
               d4 * math.log(c[2]) + d5 * math.log(c[3]))
    checks = []

    # 1) volumen
    sumL = sum(L)
    checks.append(check("determinante", abs(sumL / glv - 1) < 1e-6,
                        f"sum(L)={sumL:.3f} glv={glv:.3f} ratio={sumL/glv:.8f}"))
    # 2) GSO: pendiente del tramo inclinado vs delta_BKZ
    pendiente = -2 * math.log(delta_BKZ(b))
    tramo = [L[k + 1] - L[k] for k in range(i, min(j, len(L) - 1))]
    desv = max((abs(x - pendiente) for x in tramo), default=0.0)
    checks.append(check("GSO_perfil", desv < 0.05,
                        f"pendiente_esperada={pendiente:.6f} max_desv={desv:.6f} "
                        f"n_tramo={len(tramo)}"))
    # 3) l desde L[i]
    l = math.exp(L[i]) if L else 0.0
    checks.append(check("l_desde_GSO", l > 0 and abs(math.log(l) - L[i]) < 1e-9,
                        f"L[i]={L[i]:.4f} l={l:.2f}"))
    # 4) validez
    checks.append(check("l_menor_q", l < Q, f"l={l:.2f} q={Q}"))
    # 5..8) probabilidades y repeticiones
    deff = j - i + 1
    sigma = l / math.sqrt(deff) if deff else 0.0
    dhead1 = min(MM[0], i)
    dhead2 = min(max(0, i - dhead1), MM[1])
    dhead3 = min(max(0, i - dhead1 - dhead2), MM[2])
    dhead4 = min(max(0, i - dhead1 - dhead2 - dhead3), MM[3])
    dhead5 = min(max(0, i - dhead1 - dhead2 - dhead3 - dhead4), MM[4])
    d1 = min(MM[0] - dhead1, deff)
    d2 = min(max(0, deff - d1), MM[1] - dhead2)
    d3 = min(max(0, deff - d1 - d2), MM[2] - dhead3)
    d4 = min(max(0, deff - d1 - d2 - d3), MM[3] - dhead4)
    d5 = min(max(0, deff - d1 - d2 - d3 - d4), MM[4] - dhead5)
    ds = [d1, d2, d3, d4, d5]
    dhs = [dhead1, dhead2, dhead3, dhead4, dhead5]
    ph = [min(1.0, 2.0 * B[k] / Q) for k in range(5)]
    pm = [math.erf(c[k] * B[k + 1] / (sigma * math.sqrt(2.))) if k < 4
          else math.erf(B[4] / (sigma * math.sqrt(2.))) for k in range(5)]
    # p_middle: c aplica a B2..B5; el primero es B1 sin escalar
    pm = [math.erf(B[0] / (sigma * math.sqrt(2.)))] + [
        math.erf(c[k] * B[k + 1] / (sigma * math.sqrt(2.))) for k in range(4)]
    lg = lambda p: math.log2(p) if p > 0 else -1000.0
    log2_eps = sum(ds[k] * lg(pm[k]) + dhs[k] * lg(ph[k]) for k in range(5))
    log2_R = max(0.0, -log2_eps - nvec_sieve(b))
    svp = svp_quantum(b)
    cost = svp + log2_R
    checks.append(check("p_head", all(0 < p <= 1 for p in ph),
                        str([round(p, 6) for p in ph])))
    checks.append(check("p_middle", all(0 < p <= 1 for p in pm),
                        str([round(p, 6) for p in pm])))
    checks.append(check("epsilon", math.isfinite(log2_eps),
                        f"log2_eps={log2_eps:.3f} d={ds} dhead={dhs}"))
    checks.append(check("log2_R", log2_R >= 0,
                        f"log2_R={log2_R:.3f} nvec_sieve={nvec_sieve(b):.3f}"))
    checks.append(check("coste", math.isfinite(cost),
                        f"svp={svp:.3f} + log2_R={log2_R:.3f} = {cost:.3f}"))
    # 9) comparar contra SIS_linf_cost del estimador
    p = sis.MSISParameterSet(D, sum(M), S, *B, *M, Q, norm="linf")
    try:
        ref = sis.SIS_linf_cost(Q, w, H, B[0], MM[0], B[1], MM[1], B[2], MM[2],
                                B[3], MM[3], B[4], MM[4], b,
                                cost_svp=svp_quantum, attack_variant=v)
        dif = abs(ref - cost)
        checks.append(check("coste_vs_estimador", dif < 1.0,
                            f"reconstruido={cost:.3f} estimador={float(ref):.3f} dif={dif:.4f}"))
    except Exception as exc:
        checks.append(check("coste_vs_estimador", False, str(exc)))
    return {"variante": v, "w": w, "b": b, "h": H, "bkv": delta_BKZ(b),
            "l": round(l, 2), "sigma": round(sigma, 2),
            "coste_reconstruido": round(cost, 3),
            "checks": checks, "todo_ok": all(c["ok"] for c in checks)}


def main():
    out = {"instancia": "MatRiCT+ ASIS_est.sage", "q": Q, "d": D,
           "h": H, "variantes": {}}
    for v, (w, b) in OPTIMOS.items():
        out["variantes"][f"variant_{v}"] = validar(v, w, b)
    dest = ROOT / "security" / "validar_modelo.json"
    dest.write_text(json.dumps(out, indent=2), encoding="utf-8")
    for k, val in out["variantes"].items():
        print(f"{k}: todo_ok={val['todo_ok']} coste={val['coste_reconstruido']}")
        for c in val["checks"]:
            print(f"   [{'OK' if c['ok'] else 'FAIL'}] {c['check']}: {c['detalle']}")
    return 0 if all(v["todo_ok"] for v in out["variantes"].values()) else 1


if __name__ == "__main__":
    sys.exit(main())
