#!/usr/bin/env python3
"""Verificacion geometrica de los ataques SIS del paper (ASIS_est.sage).

Para cada variante, usando el (b,w) optimo ya encontrado, muestra:
  - volumen objetivo (det) vs volumen del perfil GSO construido
  - dimension efectiva, longitud del vector usado, condicion L<q
  - c12..c15 y como entran en el determinante
Sin optimizar: solo inspecciona la geometria de los puntos ya computados.
"""
import json
import math
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
ORIG = ROOT / "main" / "repos_originales" / "matrict_plus" / "asymmetric_sis_estimate"
sys.path.insert(0, str(ORIG))

from model_BKZ import construct_BKZ_shape, construct_BKZ_shape_randomized, delta_BKZ
import MSIS_security as sis

Q, D, SISRANK = 3489661057, 256, 5
N, LWERANK = 11, 5
B = [math.ceil(2 ** 30.3), math.ceil(2 ** 27.56), math.ceil(2 ** 22.64),
     math.ceil(2 ** 19.92), math.ceil(2 ** 16.6)]
M = [1, N + 1, SISRANK + LWERANK, 1, N + 1]
H = D * SISRANK
MM = [D * m for m in M]

# (w, b) optimos del paper_repro.json para qhat
OPTIMOS = {0: (2331, 315), 1: (2391, 350), 2: (2456, 340)}

salida = {"q": Q, "d": D, "h": H, "cotas_log2": [round(math.log2(x), 2) for x in B],
          "m_coeficientes": MM, "variantes": {}}

for variante, (w, b) in OPTIMOS.items():
    nq, n1 = w - H, H
    if variante == 0:
        c = [1.0, 1.0, 1.0, 1.0]
        i, j, L = construct_BKZ_shape(Q, nq, n1, b)
        glv = nq * math.log(Q)
    else:
        c = [B[0] / B[1], B[0] / B[2], B[0] / B[3], B[0] / B[4]] if variante == 2 else [1.0] * 4
        i, j, L = construct_BKZ_shape_randomized(Q, nq, n1, b, MM[0], MM[1], MM[2],
                                                 MM[3], MM[4], *c)
        glv = (nq * math.log(Q) + MM[1] * math.log(c[0]) + MM[2] * math.log(c[1]) +
               MM[3] * math.log(c[2]) + MM[4] * math.log(c[3]))

    lv = sum(L)
    d_ef = j - i + 1
    l = math.exp(L[i]) if L else 0.0
    sigma = l / math.sqrt(d_ef) if d_ef else 0.0

    fila = {
        "b": b, "w": w, "nq_vectores_q": nq, "n1": n1,
        "glv_log_det_objetivo": round(glv, 4),
        "sum_L_log_det_perfil": round(lv, 4),
        "ratio_det": round(lv / glv, 8) if glv else None,
        "i": i, "j": j, "dim_efectiva": d_ef,
        "c12_c13_c14_c15": [round(x, 4) for x in c],
        "rhf": round(float(delta_BKZ(b)), 6),
        "l_vec_usado": round(l, 2),
        "l_menor_q": int(l < Q),
        "sigma": round(sigma, 2),
        "coste_svp": round(sis.svp_quantum(b), 2),
    }
    salida["variantes"][f"variant_{variante}"] = fila

p = ROOT / "security" / "geom_check.json"
p.write_text(json.dumps(salida, indent=2), encoding="utf-8")
print(json.dumps(salida, indent=2))
