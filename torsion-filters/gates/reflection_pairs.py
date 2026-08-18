# -*- coding: utf-8 -*-
"""LA REFLEXION EMPAREJA LAS CLASES DE EXCESO -- verificacion de la Proposicion probada el 14.

Lo que se comprueba no es una conjetura: es una proposicion DEMOSTRADA, y este guion existe para
que la demostracion tenga su control.  El enunciado, con sigma(x) = tau - x y j* = tau - j (mod t):

  (i)   j es clase de exceso  <=>  j* lo es
  (ii)  n_j = n_{j*}
  (iii) las dos clases empatadas son exactamente los puntos fijos j = j*
  (iv)  en consecuencia e = |E| es PAR
  (v)   tau - (clase j \ {g_j}) = clase j* \ {g_{j*}}  para j no empatada

La prueba: por el Corolario 8.32, S \ g_com es estable bajo sigma, y sigma lleva la clase j a la
clase j*, luego lleva (S\g_com) inter (clase j) sobre (S\g_com) inter (clase j*).  Para j no
empatada ese conjunto es clase j \ {g_j}, no vacio, asi que j* es de exceso y los tamanos coinciden.
Las empatadas son los puntos fijos por el Corolario 8.31.  Las no fijas se emparejan, luego e es par.

CONTROL.  Las tres propiedades se comprueban tambien sobre las formas que NO se anulan.  Si se
cumplieran siempre no dirian nada de la anulacion; se cumplen en 8789 de 27235, o sea llevan
informacion.

LO QUE ESTA PROPOSICION NO DA, y es lo unico que falta para cerrar el criterio:
        tau - g_j = g_{j*}
es decir que la reflexion intercambie tambien los elementos OMITIDOS.  Con eso, tau - S = S.

Authors: Carles Marin, Claude (AI assistant).
"""
import sys
sys.path.insert(0, r"E:\proyectos\Curiosity\research\orbit-pair\gates")
from fractions import Fraction as F
from criterion_control import betas, lam_of, value, excess

PTS = [[F(3, 2), F(5, 3), F(7, 4)], [F(5, 2), F(7, 3), F(9, 4)],
       [F(4, 3), F(9, 5), F(11, 6)], [F(7, 5), F(11, 4), F(13, 7)]]


def is_zero(b, t, r):
    lam = lam_of(list(b))
    for p in PTS:
        if value(lam, t, r, p[:r]) != 0:
            return False
    return True


def data(b, t, r):
    cl = {}
    for v in b:
        cl.setdefault(v % t, []).append(v)
    if len(cl) < t:
        return None
    E = {k: sorted(v, reverse=True) for k, v in cl.items() if len(v) >= 2}
    if not E:
        return None
    S = sorted(x for c in E.values() for x in c)
    C = S[0] + S[-1]
    inc = []
    for k, c in E.items():
        inc += [(c[i] + c[i + 1], k) for i in range(len(c) - 1)]
    inc.sort(key=lambda x: -x[0])
    tau = inc[r - 1][0]
    return E, S, C, tau


CFG = [(4, 2, 16), (6, 2, 16), (8, 2, 16), (10, 2, 17), (4, 3, 14), (6, 3, 14), (2, 2, 15), (2, 3, 14)]
print("=" * 96)
print("PREDICCIONES sobre las formas que se anulan, escritas ANTES de correr:")
print("  P1  e es PAR                                  (la reflexion empareja las clases de exceso)")
print("  P2  n_j = n_{j*} para j* = tau - j (mod t)     (la biyeccion conserva el tamano)")
print("  P3  la clase j* es de exceso siempre que j lo sea")
print("  CONTROL: los mismos tres sobre las formas que NO se anulan, donde deben fallar")
print("=" * 96)
print("  t  r    Phi=0   P1 e par   P2 n_j=n_j*   P3 j* de exceso | e observados | CONTROL P1 falla")
tp1 = tp2 = tp3 = tot = 0
cbad = cn = 0
for (t, r, W) in CFG:
    nz = p1 = p2 = p3 = 0
    es = {}
    c1 = c0 = 0
    for b in betas(t, r, W):
        d = data(b, t, r)
        if d is None:
            continue
        E, S, C, tau = d
        z = is_zero(b, t, r)
        ok1 = (len(E) % 2 == 0)
        ok2 = all((tau - k) % t in E and len(E[(tau - k) % t]) == len(v) for k, v in E.items())
        ok3 = all((tau - k) % t in E for k in E)
        if not z:
            c0 += 1
            c1 += (ok1 and ok2 and ok3)
            continue
        nz += 1
        es[len(E)] = es.get(len(E), 0) + 1
        p1 += ok1
        p2 += ok2
        p3 += ok3
    tp1 += p1; tp2 += p2; tp3 += p3; tot += nz
    cbad += c1; cn += c0
    print("  %2d %2d %8d %10d %13d %17d | %-12s | %d/%d"
          % (t, r, nz, p1, p2, p3, sorted(es.items()), c1, c0))
print()
print("  P1 %d/%d    P2 %d/%d    P3 %d/%d" % (tp1, tot, tp2, tot, tp3, tot))
print("  CONTROL: las tres se cumplen tambien en %d de %d formas NO nulas" % (cbad, cn))
print("           (si fuera cn/cn, las predicciones no dirian nada de la anulacion)")
