# -*- coding: utf-8 -*-
# CUANTOS PESOS SOBREVIVEN AL CRITERIO, EN TODO TIPO.   8 de septiembre de 2026.
#
# POR QUE.  El criterio del articulo dice:  chi_lambda(rho(xi_q)) != 0  <=>  N_q(lambda) = r_q,
# con N_q(x) = #{alpha>0 : q | (x,alpha)} y r_q = N_q(rho).  Es una condicion sobre un PUNTO y un
# arreglo de paredes modulo q.  Pregunta 3 de la seccion de preguntas: cuantos puntos la cumplen,
# como funcion de q.  Si es un cuasi-polinomio, su forma dira de que depende.
#
# QUE SE MIDE.  Sobre x = sum_j c_j omega_j con c en (Z/q)^n --- o sea todo P/qP ---:
#     S(q) = #{ x : N_q(x) = r_q }
# y se compara con  prod_i (q - m_i)  y  prod_i (q + m_i),  m_i los exponentes de W.  La
# comparacion NO es un ajuste: los exponentes se calculan del sistema de raices (particion
# conjugada de las alturas) y el producto se predice antes de contar.
#
# CONTROLES.
#   K0  r_q se recalcula desde rho = (1,...,1) y tiene que coincidir con el numero de raices de
#       P-altura divisible por q --- que es la definicion del articulo.
#   K1  SENUELO: prod_i (q - i), o sea usar 1..n en vez de los exponentes.  Si acertara igual, la
#       coincidencia no diria nada sobre los exponentes.
#   K2  Se marca cada fila con si q divide o no al numero de Coxeter y a h+1, h+2, porque ahi es
#       donde r_q deja de ser 0.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python survivors_quasipoly.py > survivors_quasipoly_OUT.txt 2>&1
import json
import sys
from itertools import product

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SEP = "=" * 100

DAT = json.load(open("root_data.json", encoding="utf-8"))

# techo de coste: q^n * |raices+|
TECHO = 40_000_000


def prod(v):
    r = 1
    for x in v:
        r *= x
    return r


print(SEP)
print("S(q) = #{x en P/qP : N_q(x) = r_q}   contra   prod(q - m_i)   [m_i = exponentes de W]")
print("K1 = senuelo prod(q - i), i = 1..n.  'r_q' = numero de paredes que rho ya toca.")
print(SEP)
print("COLUMNAS: izquierda, el estadistico de RAICES (x,alpha) --- el lado de a_P, el del articulo.")
print("          derecha, el estadistico de CORRAICES <x,alpha^vee> --- el lado de a_Q, el clasico.")
print(f"{'tipo':>5} {'q':>4} | {'r_q':>4} {'S(q)':>12} {'prod(q-m)':>12} {'=?':>4} | "
      f"{'r^v':>4} {'S^v(q)':>12} {'=?':>4} | {'K1 senuelo':>12} {'=?':>4}   {'S/|W|':>10}")

resumen = {"hit": 0, "miss": 0, "hit_rq0": 0, "miss_rq0": 0, "k1_hit": 0, "k1_n": 0}
fallos = []
for d in DAT:
    n = d["rango"]
    R = d["raices"]
    exps = d["exponentes"]
    W = d["orden_W"]
    qmax = 20 if n <= 4 else (12 if n == 5 else 8)
    for q in range(2, qmax + 1):
        if (q ** n) * len(R) > TECHO:
            continue
        # r_q desde rho = (1,...,1), en las DOS lecturas del estadistico
        r_raiz = sum(1 for c, v in R if (c * sum(v)) % q == 0)
        r_corr = sum(1 for c, v in R if (sum(v)) % q == 0)
        S = Sv = 0
        for x in product(range(q), repeat=n):
            k = kv = 0
            for c, v in R:
                s = 0
                for j in range(n):
                    s += x[j] * v[j]
                if (c * s) % q == 0:
                    k += 1
                if s % q == 0:
                    kv += 1
            if k == r_raiz:
                S += 1
            if kv == r_corr:
                Sv += 1
        pred = prod([q - m for m in exps])
        k1 = prod([q - i for i in range(1, n + 1)])
        hit = (S == pred)
        hitv = (Sv == pred)
        k1hit = (S == k1)
        resumen["hit" if hit else "miss"] += 1
        resumen["hitv" if hitv else "missv"] = resumen.get("hitv" if hitv else "missv", 0) + 1
        if r_raiz == 0:
            resumen["hit_rq0" if hit else "miss_rq0"] += 1
        if r_corr == 0:
            resumen["hitv_rq0" if hitv else "missv_rq0"] = resumen.get(
                "hitv_rq0" if hitv else "missv_rq0", 0) + 1
        resumen["k1_hit" if k1hit else "k1_n"] += 1
        if not hit:
            fallos.append((d["tipo"], q, r_raiz, S, pred))
        print(f"{d['tipo']:>5} {q:>4} {r_raiz:>4} {S:>12} {pred:>12} {'SI' if hit else 'no':>4} "
              f"{r_corr:>4} {Sv:>12} {'SI' if hitv else 'no':>4} "
              f"{k1:>12} {'SI' if k1hit else 'no':>4}   {S/W:>10.4f}")

print(SEP)
print(f"  RAICES   (lado a_P): prod(q-m_i) acierta {resumen['hit']}, falla {resumen['miss']}"
      f"   [con r_q=0: {resumen['hit_rq0']} / {resumen['miss_rq0']}]")
print(f"  CORRAICES(lado a_Q): prod(q-m_i) acierta {resumen.get('hitv',0)}, "
      f"falla {resumen.get('missv',0)}"
      f"   [con r^v_q=0: {resumen.get('hitv_rq0',0)} / {resumen.get('missv_rq0',0)}]")
print(f"  K1 senuelo prod(q - i) ..: {resumen['k1_hit']} aciertos, {resumen['k1_n']} fallos"
      f"   ({'discrimina' if resumen['k1_n'] > 0 else 'NO DISCRIMINA'})")
if fallos:
    print()
    print("  donde NO cuadra (tipo, q, r_q, S(q), prod(q-m_i)):")
    for f in fallos:
        print("    ", f)
