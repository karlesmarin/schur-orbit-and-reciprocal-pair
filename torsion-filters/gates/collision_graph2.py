# -*- coding: utf-8 -*-
# EL MONOMIO TESTIGO: anatomia de las fibras AISLADAS.  15 de agosto de 2026.
#
# DE DONDE VIENE.  collision_graph.py contesta las dos primeras preguntas de la consulta externa:
#
#   1. ¿pareja unica en las formas que se anulan?   NO -- el K_2 es el 62 % (t=4 r=2), el 32 % (r=3),
#      y hay fibras de hasta 12 y 72 atomos.  La cancelacion NO tiene partner unico.
#   3. ¿vertice aislado en la familia testigo?      SI -- 80 a 120 fibras de UN SOLO atomo en cada
#      miembro que no se anula, y ninguna en los que se anulan.
#
# Un vertice aislado es un coeficiente que no puede cancelarse con nada: DEMUESTRA Phi != 0 sin mirar
# el resto de la expansion.  O sea que la ruta que propone -- "construir un monomio testigo" -- no hay
# que inventarla: los testigos ESTAN, y lo que falta es describirlos en cerrado.  Esto los describe.
#
# LO QUE SE MIDE
#   N1  la familia entera con paso 1: cuantas fibras aisladas, y sus exponentes ordenados.
#   N2  ¿es el monomio LEX-MAXIMO de la expansion uno de los aislados?  Si lo fuera, la demostracion
#       de Phi != 0 seria de una linea con un orden monomial sesgado -- que es exactamente la segunda
#       filtracion que propone la consulta.  Se prueba con los 2^r ordenes de signo y con el grado
#       total como senuelo.
#   N3  el ATOMO que produce cada testigo: su T y su reparto S.  Si el mismo (T, S) genera el testigo
#       en todos los miembros de la familia, hay formula cerrada.
#
# CONTROL
#   C1  DECOY, y arregla un confundido de collision_graph.py: alli se comparaba "% de fibras que son
#       K_2" entre las dos poblaciones, pero en una forma que se anula TODAS las fibras estan
#       equilibradas y por tanto son de tamaño PAR, mientras que en una que no se anula abundan las de
#       tamaño 1 y 3.  Parte de la diferencia era eso, no estructura.  Aqui se compara solo sobre
#       fibras EQUILIBRADAS, que es la comparacion que sí es de igual a igual.
#   C2  los miembros que se anulan tienen que dar CERO aislados.  Si dieran alguno, el instrumento
#       estaria mintiendo: un aislado es una prueba de no anulacion.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python collision_graph2.py

from collections import Counter, defaultdict

from second_stratum import setup
from peel_zero import phi_zero, C_and_tau
from collision_graph import atoms, fibras, perfil, S_of
from peel_zero import betas, occupied

t, r = 4, 2
SEM = (18, 17, 11, 8, 7, 6, 1, 0)
S0 = S_of(SEM, t)
HI, LO = S0[-1], S0[0]


def miembro(n):
    return tuple(sorted([(x + n if x == HI else (x - n if x == LO else x)) for x in SEM],
                        reverse=True))


def gcom(b):
    """los elementos de S sin reflejo en S respecto de C."""
    S = S_of(b, t)
    ct = C_and_tau(b, t, r)
    if S is None or ct is None:
        return None, None, None
    C = ct[0]
    s = set(S)
    return C, S, sorted(v for v in S if (C - v) not in s)


# GUARD __main__ -- tercera vez hoy que el mismo defecto muerde: sin el, cualquiera que
# importe una funcion de aqui se come la corrida entera por delante de su propia salida.
if __name__ == "__main__":
    print("=" * 108)
    print("N1  LA FAMILIA CON PASO 1 -- se anula y no se anula, alternando; y los testigos")
    print("=" * 108)
    print("")
    print("    n | Phi==0 |  C | g_com      | fibras | aisladas | exponentes aislados (ordenados)")
    DATA = {}
    for n in range(0, 11):
        b = miembro(n)
        if len(set(b)) != len(b):
            continue
        z = phi_zero(b, t, r)
        C, S, gc = gcom(b)
        at = atoms(b, t, r)
        if at is None:
            print("   %2d | %-6s | -- sin clases: no aplica" % (n, z))
            continue
        f = fibras(at)
        sueltas = sorted(m for m, lst in f.items() if len(lst) == 1)
        DATA[n] = (b, z, C, gc, f, sueltas)
        muestra = ', '.join(str(m) for m in sueltas[:4])
        print("   %2d | %-6s | %2d | %-10s | %6d | %8d | %s"
              % (n, z, C, str(gc), len(f), len(sueltas), muestra + (" ..." if len(sueltas) > 4 else "")))

    print("")
    print("  C2  los miembros que SE ANULAN dan cero aislados: %s"
          % ("SI, como debe ser" if all(not s for n, (_, z, _, _, _, s) in DATA.items() if z)
             else "*** NO: el instrumento miente, un aislado prueba Phi != 0 ***"))
    print("  y su g_com es %s en los que se anulan"
          % sorted({tuple(gc) for n, (_, z, _, gc, _, _) in DATA.items() if z}))
    print("  frente a %s en los que no"
          % sorted({tuple(gc) for n, (_, z, _, gc, _, _) in DATA.items() if not z}))

    # ===================================================================== N2 ========================
    print("")
    print("=" * 108)
    print("N2  ¿ES EL EXTREMO DE ALGUN ORDEN MONOMIAL UN TESTIGO AISLADO?")
    print("=" * 108)
    print("")
    print("  Si lo es, 'Phi != 0' sale de mirar UN termino, que es la segunda filtracion que se propone.")
    print("")
    print("    n | grado total (señuelo)      | lex z1>>z2      | lex z1>>-z2     | -z1>>z2")
    for n in sorted(DATA):
        b, z, C, gc, f, sueltas = DATA[n]
        if z:
            continue
        ss = set(sueltas)
        filas = []
        # senuelo: el grado total maximo -- ya sabemos que su estrato se anula, luego debe FALLAR
        dmax = max(sum(m) for m in f)
        cand = [m for m in f if sum(m) == dmax]
        filas.append("%-26s" % ("%d cand, aislado? %s" % (len(cand), any(m in ss for m in cand))))
        for signos in [(1, 1), (1, -1), (-1, 1)]:
            key = lambda m: tuple(s * x for s, x in zip(signos, m))
            top = max(f, key=key)
            filas.append("%-15s" % ("%s %s" % (top, "AISLADO" if top in ss else "no")))
        print("   %2d | %s | %s | %s | %s" % (n, filas[0], filas[1], filas[2], filas[3]))

    # ===================================================================== N3 ========================
    print("")
    print("=" * 108)
    print("N3  EL ATOMO QUE PRODUCE EL TESTIGO -- ¿es el mismo (T, S) en toda la familia?")
    print("=" * 108)
    print("")
    print("    n | testigo lex-minimo | T del atomo                    | reparto S | rel. con C y g_com")
    for n in sorted(DATA):
        b, z, C, gc, f, sueltas = DATA[n]
        if z or not sueltas:
            continue
        m = sueltas[0]
        (c, tag) = f[m][0]
        T, Sp = tag
        print("   %2d | %-18s | %-30s | %-9s | m1+C=%d  m2=%s"
              % (n, str(m), str(T), str(Sp), m[0] + C, m[1]))

    print("")
    print("=" * 108)
    print("C1  DECOY CORREGIDO -- solo fibras EQUILIBRADAS, que es la comparacion de igual a igual")
    print("=" * 108)
    print("")
    print("   t  r   W | poblacion            n | fibras equilibradas | de ellas, de tamaño 2")
    for (tt, rr, W) in [(4, 2, 14), (6, 2, 15)]:
        for quiero_cero in (True, False):
            n = eq = eq2 = 0
            for b in betas(tt, rr, W):
                if not occupied(b, tt):
                    continue
                zz = phi_zero(b, tt, rr)
                if zz is None or bool(zz) != quiero_cero:
                    continue
                at = atoms(b, tt, rr)
                if at is None:
                    continue
                n += 1
                if n > 40:
                    break
                for m, lst in fibras(at).items():
                    p = sum(1 for c, _ in lst if c > 0)
                    if p * 2 == len(lst):
                        eq += 1
                        eq2 += (len(lst) == 2)
            et = "Phi == 0" if quiero_cero else "Phi != 0  [DECOY]"
            if not n:
                print("  %3d %2d %3d | %-21s %2d | POBLACION VACIA" % (tt, rr, W, et, 0))
                continue
            print("  %3d %2d %3d | %-21s %2d | %19d | %5.1f %%"
                  % (tt, rr, W, et, n, eq, 100.0 * eq2 / max(eq, 1)))
    print("")
    print("  Si los dos porcentajes se parecen, el K_2 del guion anterior era en buena parte el efecto")
    print("  de que en una forma que se anula NO HAY fibras impares, y no una diferencia de estructura.")
    print("")
    print("=" * 108)
