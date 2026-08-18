# -*- coding: utf-8 -*-
# LOS VERTICES DEL POLIGONO DE NEWTON SON ATOMOS SUELTOS.  15 de agosto de 2026.
#
# DE DONDE SALE, y es una correccion a mi propia medida.  collision_graph2.py pregunto si el monomio
# extremo de algun orden monomial es un testigo aislado y contesto QUE NO en los tres ordenes
# probados.  Ese NO estaba MAL MEDIDO: tomaba el extremo sobre TODAS las fibras, incluidas las que se
# cancelan.  El extremo de Phi no es el extremo de los candidatos: es el extremo de lo que SOBREVIVE.
# Rehecho sobre el soporte real, el resultado se da la vuelta entera.
#
# LO QUE SE MIDE
#   N1  la familia testigo: el poligono de Newton de Phi y cuantos atomos tiene cada vertice.
#   N2  la POBLACION: en las formas con Phi != 0, ¿es TODO vertice del poligono de Newton un atomo
#       suelto?  Si lo es en general, "Phi != 0" se certifica exhibiendo un vertice, y la segunda
#       filtracion que propone la consulta externa -- un orden monomial sesgado -- existe de verdad.
#   N3  la forma del poligono contra invariantes de S, para ver si los vertices se describen cerrados.
#
# EL TEST EN CUALQUIER DIMENSION.  Para r > 2 no se calcula la envolvente: se maximizan funcionales
# lineales enteros al azar sobre el soporte.  El argmax de un funcional generico ES un vertice, y si
# el argmax tiene un solo atomo, ese vertice es un testigo.  Sirve para r = 2 y para r = 3 igual.
#
# CONTROLES, y los dos pueden fallar -- de hecho uno TIENE que fallar
#   C1  DECOY OBLIGADO: los vertices de la envolvente de TODAS las fibras (las que cancelan tambien).
#       Si esos tambien fueran siempre sueltos, Phi no podria anularse NUNCA y el instrumento estaria
#       midiendo una tautologia.  Este control tiene que FALLAR, y si no falla, N2 no vale.
#   C2  poblacion: se imprime n SIEMPRE, y las formas que se anulan se excluyen diciendolo -- su
#       soporte es vacio y no tienen poligono.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python newton_vertices.py

import random
from collections import Counter

from collision_graph import atoms, fibras, S_of
from peel_zero import betas, occupied, phi_zero, C_and_tau

RNG = random.Random(20260815)


def hull2(P):
    """envolvente convexa en el plano, sentido antihorario, sin puntos colineales."""
    P = sorted(set(P))
    if len(P) < 3:
        return list(P)

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    lo = []
    for p in P:
        while len(lo) >= 2 and cross(lo[-2], lo[-1], p) <= 0:
            lo.pop()
        lo.append(p)
    up = []
    for p in reversed(P):
        while len(up) >= 2 and cross(up[-2], up[-1], p) <= 0:
            up.pop()
        up.append(p)
    return lo[:-1] + up[:-1]


def soporte(f):
    return [m for m, lst in f.items() if sum(c for c, _ in lst) != 0]


def extremos_azar(pts, r, k=40):
    """argmax de k funcionales lineales genericos: cada uno es un vertice."""
    out = set()
    for _ in range(k):
        w = tuple(RNG.randint(-1000, 1000) for _ in range(r))
        best = max(pts, key=lambda m: sum(a * b for a, b in zip(w, m)))
        out.add(best)
    return out


# ===================================================================== N1 ========================
# GUARD __main__ -- tercera vez hoy que el mismo defecto muerde: sin el, cualquiera que
# importe una funcion de aqui se come la corrida entera por delante de su propia salida.
if __name__ == "__main__":
    t, r = 4, 2
    SEM = (18, 17, 11, 8, 7, 6, 1, 0)
    S0 = S_of(SEM, t)
    HI, LO = S0[-1], S0[0]

    print("=" * 108)
    print("N1  LA FAMILIA TESTIGO -- el poligono de Newton de Phi, vertice a vertice")
    print("=" * 108)
    print("")
    print("    n | Phi==0 | soporte | vertices | atomos por vertice | los vertices")
    for n in range(0, 11):
        b = tuple(sorted([(x + n if x == HI else (x - n if x == LO else x)) for x in SEM], reverse=True))
        if len(set(b)) != len(b):
            continue
        z = phi_zero(b, t, r)
        at = atoms(b, t, r)
        if at is None:
            continue
        f = fibras(at)
        sop = soporte(f)
        if not sop:
            print("   %2d | %-6s |       0 | -- se anula: no hay poligono" % (n, z))
            continue
        V = hull2(sop)
        cuenta = sorted(len(f[v]) for v in V)
        print("   %2d | %-6s | %7d | %8d | %-18s | %s"
              % (n, z, len(sop), len(V), str(cuenta), str(sorted(V))[:44]))

    print("")
    print("  Un vertice con UN atomo es un coeficiente +-1 que no puede cancelarse con nada:")
    print("  demuestra Phi != 0 mirando un solo termino.")

    # ===================================================================== N2 + C1 ===================
    print("")
    print("=" * 108)
    print("N2  LA POBLACION -- ¿es TODO vertice un atomo suelto?   C1  y el DECOY que tiene que fallar")
    print("=" * 108)
    print("")
    print("   t  r   W |   n | N2: formas con TODOS los vertices sueltos | C1 decoy: idem con las fibras")
    print("             |     |     (soporte real de Phi)                |     que cancelan incluidas")
    for (tt, rr, W) in [(4, 2, 14), (4, 2, 15), (6, 2, 15), (4, 3, 14), (6, 3, 16)]:
        n = ok = ok_decoy = 0
        peor = None
        for b in betas(tt, rr, W):
            if not occupied(b, tt):
                continue
            z = phi_zero(b, tt, rr)
            if z is None or z:
                continue                       # las que se anulan no tienen poligono: C2
            at = atoms(b, tt, rr)
            if at is None:
                continue
            n += 1
            if n > 120:
                break
            f = fibras(at)
            sop = soporte(f)
            if not sop:
                continue
            V = hull2(sop) if rr == 2 else extremos_azar(sop, rr)
            c = [len(f[v]) for v in V]
            if all(x == 1 for x in c):
                ok += 1
            elif peor is None:
                peor = (b, sorted(c))
            Vd = hull2(list(f)) if rr == 2 else extremos_azar(list(f), rr)
            if all(len(f[v]) == 1 for v in Vd):
                ok_decoy += 1
        if not n:
            print("  %3d %2d %3d | %3d | POBLACION VACIA -- esta fila no dice nada" % (tt, rr, W, 0))
            continue
        print("  %3d %2d %3d | %3d | %5d  (%5.1f %%)                          | %5d  (%5.1f %%)"
              % (tt, rr, W, n, ok, 100.0 * ok / n, ok_decoy, 100.0 * ok_decoy / n))
        if peor:
            print("             |     | primer contraejemplo: %s  atomos %s" % (peor[0], peor[1]))

    print("")
    print("  C1  el decoy TIENE que dar un porcentaje mas bajo: si los vertices de TODAS las fibras")
    print("      fueran siempre sueltos, Phi no podria anularse jamas y N2 seria una tautologia.")

    # ===================================================================== C3 ========================
    print("")
    print("=" * 108)
    print("C3  EL CONTROL DONDE SI PUEDE FALLAR -- sobre las formas que SE ANULAN")
    print("=" * 108)
    print("")
    print("  Si Phi == 0, TODO vertice de la envolvente de los exponentes CANDIDATOS tiene que llevar")
    print("  >= 2 atomos: un vertice con uno solo seria un +-1 que nada cancela, o sea Phi != 0.  Este")
    print("  control PUEDE fallar y su fallo seria una contradiccion, no un matiz.")
    print("")
    print("   t  r   W |   n | vertices con 1 solo atomo | minimo de atomos por vertice")
    for (tt, rr, W) in [(4, 2, 14), (6, 2, 15), (4, 3, 14), (6, 3, 16)]:
        n = viol = 0
        mins = []
        for b in betas(tt, rr, W):
            if not occupied(b, tt):
                continue
            z = phi_zero(b, tt, rr)
            if z is None or not z:
                continue
            at = atoms(b, tt, rr)
            if at is None:
                continue
            n += 1
            if n > 80:
                break
            f = fibras(at)
            V = hull2(list(f)) if rr == 2 else extremos_azar(list(f), rr)
            c = [len(f[v]) for v in V]
            if not c:
                continue
            mins.append(min(c))
            viol += sum(1 for x in c if x == 1)
        if not n:
            print("  %3d %2d %3d | %3d | POBLACION VACIA -- esta fila no dice nada" % (tt, rr, W, 0))
            continue
        print("  %3d %2d %3d | %3d | %25d | %d   %s"
              % (tt, rr, W, n, viol, min(mins) if mins else -1,
                 "consistente" if viol == 0 else "*** IMPOSIBLE: revisar el instrumento ***"))
    print("")
    print("  Esto es la version CON CONTENIDO de la 'residual rigidity' que se propone: la anulacion")
    print("  no obliga a una involucion, pero SI obliga a colision en los EXTREMOS.")

    # ===================================================================== N3 ========================
    print("")
    print("=" * 108)
    print("N3  ¿SE DESCRIBEN LOS VERTICES EN CERRADO?  contra invariantes de S")
    print("=" * 108)
    print("")
    print("    n | S                                    | C  | vertices (a,b) con |a|>=|b| | a vs maxS-1")
    for n in range(0, 11, 2):
        b = tuple(sorted([(x + n if x == HI else (x - n if x == LO else x)) for x in SEM], reverse=True))
        at = atoms(b, t, r)
        f = fibras(at)
        sop = soporte(f)
        if not sop:
            continue
        V = hull2(sop)
        S = S_of(b, t)
        C = C_and_tau(b, t, r)[0]
        ab = sorted({(max(abs(x), abs(y)), min(abs(x), abs(y))) for x, y in V})
        a = ab[0][0] if len(ab) == 1 else None
        print("   %2d | %-36s | %2d | %-27s | maxS-1 = %d  %s"
              % (n, str(S)[:36], C, str(ab), max(S) - 1,
                 "IGUAL" if a == max(S) - 1 else ("distinto" if a else "varios")))

    print("")
    print("=" * 108)
