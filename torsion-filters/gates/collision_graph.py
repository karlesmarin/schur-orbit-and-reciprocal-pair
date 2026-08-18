# -*- coding: utf-8 -*-
# EL GRAFO DE COLISIONES: ¿la cancelacion de Phi tiene PAREJA UNICA?  15 de agosto de 2026.
#
# DE DONDE SALE.  La consulta externa, vuelta 06, contesta P6 con un NO honesto (no hay teorema de
# extension de involuciones que sirva) y propone en su lugar reformular la Conjetura 8.44 asi:
#
#     RESIDUAL RIGIDITY:  g_com asimetrico  =>  alguna FIBRA DE EXPONENTE esta desequilibrada
#
# El puente es elemental y correcto: Phi es un polinomio de Laurent en z, luego Phi == 0 obliga a que
# la suma de coeficientes se anule EXPONENTE A EXPONENTE.  (El invoca independencia de caracteres de
# Dedekind; en nuestro caso es mas simple todavia -- los monomios son una base -- y asi lo decimos.)
# Y propone MEDIR antes de intentar demostrar nada, con este grafo:
#
#     vertices = atomos de la expansion       arista  omega ~ omega'  si  m_omega = m_omega'  y
#                                                                        c_omega = -c_omega'
#
# y tres preguntas suyas, que son las tres columnas N1-N3:
#
#   1. en las formas que SE ANULAN, ¿es cada componente un K_2?  (=> la pareja es UNICA)
#   2. si lo es, ¿la pareja corresponde a sustituir el dato de g_com por su reflejo C-x?
#   3. en la familia testigo que NO se anula, ¿aparece una fibra desequilibrada -- idealmente un
#      VERTICE AISLADO?  Ese seria el monomio testigo que su ruta necesita construir.
#
# LOS ATOMOS SON EXACTAMENTE LOS DE stratum().  No hay reimplementacion: el mismo desarrollo de
# Laplace por bloques, el mismo split_sign, los mismos alternantes.  Aqui solo se deja de sumar y se
# guarda cada sumando por separado.  Y sus coeficientes SON +-1 exactos: w = perm_sign es +-1,
# split_sign es +-1 y los coeficientes del alternante son signos de permutacion.  Eso hace que la
# pregunta de la pareja tenga sentido literal; si fueran enteros generales, "tantos + como -" no
# seria una condicion de emparejamiento.
#
# CONTROLES, y los tres pueden fallar
#   C0  FATAL.  La suma de los atomos por exponente reproduce la expansion entera de full_dict(), y
#       "todas las fibras se anulan" coincide con phi_zero().  Si no, los atomos no son Phi.
#   C1  DECOY, y es el control que este guion existe para no saltarse.  La misma estadistica de
#       fibras sobre formas que NO se anulan.  Si el histograma sale parecido, la estadistica no
#       distingue las dos poblaciones y N1 no mide lo que dice medir.  Un empate no es un acuerdo.
#   C2  no vacuidad: se imprime SIEMPRE cuantas formas hay en cada poblacion.  Con cero formas que se
#       anulan no se dice nada de la pareja unica.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python collision_graph.py

import itertools
from collections import Counter, defaultdict
from itertools import combinations

from second_stratum import setup, all_transversals
from depth_histogram import alt, split_sign, full_dict
from peel_zero import betas, occupied, phi_zero, C_and_tau


def atoms(beta, t, r):
    """[(m, c, tag)] -- un atomo por sumando del desarrollo de Laplace, c = +-1."""
    st = setup(beta, t)
    if st is None:
        return None
    cl, E, Cd = st
    if not E:
        return None
    tr = all_transversals(beta, cl, r, t)
    idx = tuple(range(2 * r))
    out = []
    for (_, T, w, _) in tr:
        for S in combinations(idx, r):
            Sc = tuple(a for a in idx if a not in S)
            base = w * split_sign(S, Sc, r)
            A = alt([T[a] for a in S], r)
            B = alt([-T[a] for a in Sc], r)
            for ka, ca in A.items():
                for kb, cb in B.items():
                    m = tuple(ka[j] + kb[j] for j in range(r))
                    out.append((m, base * ca * cb, (T, S)))
    return out


def fibras(at):
    f = defaultdict(list)
    for (m, c, tag) in at:
        f[m].append((c, tag))
    return f


def perfil(f):
    """(n_fibras, hist de (n+, n-), n_K2, n_desequilibradas, n_aisladas)."""
    hist = Counter()
    k2 = des = ais = 0
    for m, lst in f.items():
        p = sum(1 for c, _ in lst if c > 0)
        q = len(lst) - p
        hist[(p, q)] += 1
        if (p, q) == (1, 1):
            k2 += 1
        if p != q:
            des += 1
        if len(lst) == 1:
            ais += 1
    return len(f), hist, k2, des, ais


def S_of(beta, t):
    st = setup(beta, t)
    if st is None:
        return None
    cl, E, Cd = st
    if not E:
        return None
    return sorted({v for k in E for v in Cd[k]})


# GUARD __main__ -- y esto es el MISMO defecto que este guion documenta haber arreglado en
# flag.py hoy mismo, cometido acto seguido por su propio autor: collision_graph2.py importa
# atoms() de aqui, y sin guard la corrida entera se colaba por delante de la suya.  Cambio de
# pura indentacion, verificado por diff inverso.
if __name__ == "__main__":
    # ===================================================================== C0 ========================
    print("=" * 108)
    print("C0  ACEPTACION -- los atomos SON Phi, y su balance por fibra ES la anulacion")
    print("=" * 108)
    malo = 0
    vistos = 0
    for (t, r, W) in [(4, 2, 13), (6, 2, 14)]:
        for b in betas(t, r, W):
            if not occupied(b, t):
                continue
            at = atoms(b, t, r)
            if at is None:
                continue
            vistos += 1
            if vistos > 400:
                break
            acc = defaultdict(int)
            for (m, c, _) in at:
                acc[m] += c
            mio = {k: v for k, v in acc.items() if v}
            st = setup(b, t)
            tr = all_transversals(b, st[0], r, t)
            ref = full_dict([(x[2], x[1]) for x in tr], r)
            if mio != ref:
                malo += 1
            if (not mio) != bool(phi_zero(b, t, r)):
                malo += 1
    print("     %d formas comparadas atomo a atomo contra full_dict() y contra phi_zero()" % min(vistos, 400))
    print("     C0 %s" % ("PASA" if not malo else "*** FALLA en %d ***" % malo))
    if malo:
        print("DONE (veredicto suspendido)")
        raise SystemExit(1)

    # ===================================================================== N1 + C1 ===================
    print("")
    print("=" * 108)
    print("N1  LAS FORMAS QUE SE ANULAN  --  y C1, el DECOY: las que no, con la misma estadistica")
    print("=" * 108)
    print("")
    print("   t  r   W | poblacion             n | fibras/forma | %K_2  | tam. de fibra mas comun")
    for (t, r, W) in [(4, 2, 14), (4, 2, 15), (6, 2, 15), (4, 3, 15)]:
        for quiero_cero in (True, False):
            n = 0
            tot_f = tot_k2 = 0
            tam = Counter()
            maxfib = 0
            for b in betas(t, r, W):
                if not occupied(b, t):
                    continue
                z = phi_zero(b, t, r)
                if z is None or bool(z) != quiero_cero:
                    continue
                at = atoms(b, t, r)
                if at is None:
                    continue
                n += 1
                if n > 60:
                    break
                f = fibras(at)
                nf, hist, k2, des, ais = perfil(f)
                tot_f += nf
                tot_k2 += k2
                for (p, q), c in hist.items():
                    tam[p + q] += c
                    maxfib = max(maxfib, p + q)
            etiqueta = "Phi == 0 (conj 8.44)" if quiero_cero else "Phi != 0  [DECOY C1]"
            if not n:
                print("  %3d %2d %3d | %-20s %3d | -- POBLACION VACIA: no se dice nada de esta fila"
                      % (t, r, W, etiqueta, 0))
                continue
            comunes = ', '.join('%d:%d' % (k, v) for k, v in tam.most_common(4))
            print("  %3d %2d %3d | %-20s %3d | %12.1f | %5.1f | %s   (max %d)"
                  % (t, r, W, etiqueta, n, tot_f / float(n), 100.0 * tot_k2 / max(tot_f, 1),
                     comunes, maxfib))

    print("")
    print("  LECTURA.  Si el %K_2 de la fila Phi==0 no es 100, la cancelacion NO tiene pareja unica y su")
    print("  pregunta 1 se contesta que NO.  Si ademas el DECOY da un perfil parecido, la estadistica no")
    print("  separa las dos poblaciones y no hay nada que leer en la primera fila.")

    # ===================================================================== N2 ========================
    print("")
    print("=" * 108)
    print("N2  SU PREGUNTA 3 -- la familia testigo que NO se anula: ¿hay fibra desequilibrada AISLADA?")
    print("=" * 108)
    SEMILLA = (18, 17, 11, 8, 7, 6, 1, 0)
    t, r = 4, 2
    S0 = S_of(SEMILLA, t)
    hi, lo = S0[-1], S0[0]
    print("")
    print("  familia:  max S sube n, min S baja n,  semilla %s" % (SEMILLA,))
    print("")
    print("     n |  W | fibras | desequilibradas | aisladas | exponentes de las aisladas (hasta 3)")
    for n in range(0, 9):
        b = tuple(sorted([(x + n if x == hi else (x - n if x == lo else x)) for x in SEMILLA],
                         reverse=True))
        if len(set(b)) != len(b):
            continue
        at = atoms(b, t, r)
        if at is None:
            continue
        f = fibras(at)
        nf, hist, k2, des, ais = perfil(f)
        sueltas = [m for m, lst in f.items() if len(lst) == 1]
        muestra = ', '.join(str(m) for m in sorted(sueltas)[:3]) or "(ninguna)"
        print("   %3d | %2d | %6d | %15d | %8d | %s" % (n, b[0] - b[-1], nf, des, ais, muestra))

    print("")
    print("  Si la columna 'aisladas' es > 0 y las mismas fibras siguen solas al crecer n, el monomio")
    print("  testigo EXISTE en la familia y se puede intentar describir en cerrado.  Si es 0, su ruta")
    print("  necesita un testigo que no sea un vertice aislado, y eso es mas trabajo, no menos.")
    print("")
    print("=" * 108)
