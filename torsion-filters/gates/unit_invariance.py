# -*- coding: utf-8 -*-
# EL LUGAR DEL FILTRO ES GALOIS-INVARIANTE.   16 de agosto de 2026.
#
# DE DONDE SALE.  filter_translate.py encontro que a orden impar el lugar regular de tipo B es el de
# tipo C trasladado por s = 2^{-1}.  Reescrito en exponentes, la traslacion NO es una traslacion:
#
#     b = 2u+1  (tipo B),   c = u+s  =>  c = s(2u+1) = s.b .
#
# O sea, el vector de exponentes de tipo C es el de tipo B MULTIPLICADO POR LA UNIDAD 2^{-1}.
# Y multiplicar por una unidad de Z/t es un automorfismo que preserva las tres condiciones de
# regularidad (x != 0, x_i != x_j, x_i != -x_j), porque las preserva TODAS a la vez.
#
# LA CONJETURA GRANDE, que es la que hay que medir:
#
#   (U)  el lugar regular en coordenadas de exponente es invariante bajo TODO el grupo de unidades
#        (Z/t)^x, no solo bajo 2^{-1}.
#
# Y si (U) vale, la lectura no es aritmetica sino de GALOIS.  Gal(Q(zeta_t)/Q) = (Z/t)^x actua por
# zeta -> zeta^k; el filtro tau toma valores en {0,+-1} c Q, luego es fijo por Galois, luego su
# lugar de anulacion tiene que ser (Z/t)^x-invariante.  La invariancia no es una coincidencia: es la
# RACIONALIDAD del filtro, que el paper ya usa sin nombrarla.
#
# Y entonces la "segunda manifestacion de la paridad" se lee de golpe:
#
#     los filtros de tipo B y de tipo C ven el mismo lugar  <=>  2 es una unidad  <=>  t impar,
#
# porque el paso de uno al otro es multiplicar por 2, que es Galois exactamente cuando es unidad.
#
# LO QUE SE MIDE
#   U1  para cada t y cada m: ¿es el lugar regular invariante bajo cada k en (Z/t)^x?
#   U2  y bajo los k NO invertibles (los que existen si t no es primo): tienen que FALLAR, o la
#       invariancia no dice nada sobre las unidades en particular.
#   U3  el signo: signo(v) contra signo(k.v).  ¿es constante en v para cada k? -- si lo fuera, el
#       filtro entero (no solo su lugar) seria equivariante.
#   U4  EL RENDIMIENTO.  Si (U) vale, el filtro no hay que evaluarlo en todo el lugar sino en UN
#       PUNTO POR ORBITA de (Z/t)^x.  Se cuentan las orbitas y se da el factor |lugar|/|orbitas|,
#       que es lo que se ahorra de verdad -- y se compara con phi(t), que es la cota de arriba
#       (las orbitas libres ahorran phi(t); las que tienen estabilizador, menos).
#
# CONTROLES
#   C0  U2 es el control que puede fallar: si los NO invertibles tambien preservaran el lugar, la
#       propiedad seria trivial y no distinguiria nada.
#   C1  se cuenta el lugar bruto y se imprime; los casos vacios se marcan VACIO y no cuentan.
#   C2  n impreso siempre.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python unit_invariance.py
import itertools
import json
from math import gcd


def regular(v, t):
    cl = []
    for x in v:
        cl.append(x % t)
        cl.append((-x) % t)
    return (0 not in cl) and len(set(cl)) == len(cl)


def signo(v, t):
    val, sg = [], 1
    for x in v:
        e = x % t
        if e > t - e:
            e, s = t - e, -1
        else:
            s = 1
        val.append(e)
        sg *= s
    idx = sorted(range(len(val)), key=lambda i: -val[i])
    vis = [False] * len(idx)
    for i in range(len(idx)):
        if vis[i]:
            continue
        j, ciclo = i, 0
        while not vis[j]:
            vis[j] = True
            j = idx[j]
            ciclo += 1
        if ciclo % 2 == 0:
            sg = -sg
    return sg


print("=" * 118)
print("(U)  ¿ES EL LUGAR DEL FILTRO INVARIANTE BAJO EL GRUPO DE UNIDADES?")
print("=" * 118)

RES = []
tot_u = ok_u = 0
tot_nu = ok_nu = 0
for t in range(3, 13):
    for m in [1, 2, 3]:
        if t ** m > 200000:
            continue
        TODOS = list(itertools.product(range(t), repeat=m))
        reg = set(v for v in TODOS if regular(v, t))
        if not reg:
            print("  t=%2d m=%d : lugar VACIO, no mide" % (t, m))
            continue
        unidades, nounidades = [], []
        signos_ctes = []
        for k in range(1, t):
            img = set(tuple((k * x) % t for x in v) for v in reg)
            inv = (img == reg)
            if gcd(k, t) == 1:
                unidades.append((k, inv))
                tot_u += 1
                ok_u += inv
                if inv:
                    q = set(signo(v, t) * signo([(k * x) % t for x in v], t) for v in reg)
                    signos_ctes.append((k, sorted(q)))
            else:
                nounidades.append((k, inv))
                tot_nu += 1
                ok_nu += inv
        malas_u = [k for k, i in unidades if not i]
        buenas_nu = [k for k, i in nounidades if i]
        cte = all(len(q) == 1 for _, q in signos_ctes)
        # U4 orbitas del grupo de unidades sobre el lugar
        U = [k for k in range(1, t) if gcd(k, t) == 1]
        vistos, orbitas = set(), 0
        for v in sorted(reg):
            if v in vistos:
                continue
            orbitas += 1
            for k in U:
                vistos.add(tuple((k * x) % t for x in v))
        phi = len(U)
        print("")
        print("  t=%2d m=%d   |lugar| = %d de %d" % (t, m, len(reg), len(TODOS)))
        print("     U1  unidades que PRESERVAN el lugar : %d de %d   %s"
              % (len(unidades) - len(malas_u), len(unidades),
                 "todas" if not malas_u else "!! fallan %s" % malas_u))
        print("     U2  NO-unidades que preservan el lugar : %s  (%d en total)"
              % (buenas_nu if buenas_nu else "ninguna", len(nounidades)))
        print("     U3  signo(v).signo(k.v) constante en v, para cada unidad : %s"
              % ("si, valores %s" % [q for _, q in signos_ctes] if cte else "NO"))
        print("     U4  orbitas de Galois sobre el lugar : %d  ->  factor %.2f  (phi(t) = %d)"
              % (orbitas, len(reg) / float(orbitas), phi))
        RES.append({"t": t, "m": m, "n_lugar": len(reg), "n_total": len(TODOS),
                    "unidades_que_fallan": malas_u,
                    "no_unidades_que_preservan": buenas_nu,
                    "n_no_unidades": len(nounidades),
                    "signo_constante": bool(cte),
                    "n_orbitas": int(orbitas), "phi": int(phi),
                    "factor": round(len(reg) / float(orbitas), 3)})

print("")
print("=" * 118)
print("  RESUMEN   unidades que preservan : %d de %d" % (ok_u, tot_u))
print("            NO-unidades que preservan : %d de %d   <- el control; tiene que ser bajo"
      % (ok_nu, tot_nu))
print("  LECTURA, escrita ANTES de correr:")
print("   * si TODAS las unidades preservan el lugar y casi ninguna no-unidad lo hace, (U) vale y")
print("     el lugar del filtro es Galois-invariante -- y B/C coinciden a orden impar porque 2 es")
print("     unidad, no por una coincidencia de paredes.")
print("   * si alguna unidad fallara, (U) es falsa y la coincidencia B/C es solo con el 2.")
json.dump(RES, open("unit_invariance_DUMP.json", "w"), indent=1)
print("=" * 118)
print("DONE")
