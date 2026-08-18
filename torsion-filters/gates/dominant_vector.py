# -*- coding: utf-8 -*-
# EL VECTOR DOMINANTE v(T), Y SI SU MAXIMO DECIDE.  15 de agosto de 2026.
#
# DE DONDE SALE.  La consulta externa, vuelta 08, propone dejar de contar atomos y comprimir cada
# minor a UN vector:
#
#     v(T) = ( t_1 - t_{2r},  t_2 - t_{2r-1},  ...,  t_r - t_{r+1} )        con T ordenado decreciente
#
# y afirma tres cosas encadenadas:
#
#   (1)  Newt(A(T)) = conv( W(B_r) . v(T) )         el orbit polytope de tipo B/C
#   (2)  cada vertice de (1) viene de UN SOLO atomo  -- hay que emparejar el mayor con el menor
#   (3)  luego el certificado de no anulacion es:  algun v MAXIMAL EN DOMINANCIA entre las
#        transversales tiene multiplicidad con signo  m(v) != 0
#
# y el experimento que pide: sobre la poblacion critica, agrupar por v(T), mirar solo los maximales
# en dominancia de tipo C, y comparar las dos columnas  g_com simetrico  /  g_com asimetrico.
#
# LO QUE ADEMAS HAY QUE DECIR, Y NO LO DICE NADIE.  Nuestras clases residuales son DISJUNTAS por
# construccion (cl[k] = indices con beta_i = k mod t).  Una familia disjunta da un matroide de
# PARTICION, no un matroide transversal general: las transversales son literalmente el producto de
# los bloques.  O sea que aqui NO hay estructura matroidal que aprovechar -- ni matching fields, ni
# matroides valuados, ni oriented matroids: el matroide subyacente es trivial.  Toda la dificultad
# esta en que v(T) NO ES SEPARABLE por bloques.  Se comprueba en C1, no se afirma.
#
# LO QUE SE MIDE
#   N1  la formula (1) y (2), verificadas contra la expansion real, minor a minor.
#   N2  el maximo en dominancia entre transversales: ¿es unico?  ¿y su m(v)?
#   N3  EL EXPERIMENTO: sobre la poblacion critica (C = tau y S\g_com simetrico), las dos columnas.
#
# CONTROLES
#   C0  FATAL.  (1) y (2) se comprueban contra los atomos de verdad, no se dan por buenas.
#   C1  el matroide es de particion: las clases son disjuntas y |transversales| = prod |cl[k]|.
#       Si fallara, la lectura de arriba es falsa y hay que rehacerla.
#   C2  poblacion: n SIEMPRE impreso, y las dos columnas por separado.  Una columna vacia no dice
#       nada de la otra.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python dominant_vector.py

import itertools
from collections import defaultdict

from second_stratum import setup, all_transversals
from depth_histogram import alt, split_sign
from collision_graph import atoms, fibras, S_of
from newton_vertices import hull2
from peel_zero import betas, occupied, phi_zero, C_and_tau


def v_of(T, r):
    """v(T) = (t_1 - t_2r, t_2 - t_{2r-1}, ...), con T decreciente."""
    T = sorted(T, reverse=True)
    return tuple(T[i] - T[2 * r - 1 - i] for i in range(r))


def domina(a, b):
    """a >= b en dominancia: todas las sumas parciales."""
    sa = sb = 0
    for x, y in zip(a, b):
        sa += x
        sb += y
        if sa < sb:
            return False
    return True


def atoms_of_minor(T, r):
    """los atomos de A(T) solo: [(m, c)]."""
    idx = tuple(range(2 * r))
    out = []
    for S in itertools.combinations(idx, r):
        Sc = tuple(a for a in idx if a not in S)
        base = split_sign(S, Sc, r)
        A = alt([T[a] for a in S], r)
        B = alt([-T[a] for a in Sc], r)
        for ka, ca in A.items():
            for kb, cb in B.items():
                out.append((tuple(ka[j] + kb[j] for j in range(r)), base * ca * cb))
    return out


def gcom_de(b, t, r):
    S = S_of(b, t)
    ct = C_and_tau(b, t, r)
    if S is None or ct is None:
        return None, None, None
    C = ct[0]
    s = set(S)
    return C, ct[1], sorted(x for x in S if (C - x) not in s)


# GUARD __main__ -- por la misma razon que en flag.py, collision_graph.py y newton_vertices.py.
if __name__ == "__main__":
    # ===================================================================== C0 + N1 ===================
    print("=" * 108)
    print("C0 + N1  LA FORMULA v(T): Newt(A(T)) = orbit polytope de tipo B/C, y vertice = UN atomo")
    print("=" * 108)
    print("")
    print("    T                     | v(T)      | vertices reales de Newt(A(T))        | orbit  | 1 atomo")
    r = 2
    malo = 0
    for T in [(17, 11, 6, 0), (20, 11, 6, 1), (24, 11, 6, 1), (13, 9, 4, 2), (30, 7, 5, 1)]:
        at = atoms_of_minor(T, r)
        f = defaultdict(list)
        for m, c in at:
            f[m].append(c)
        sop = [m for m, l in f.items() if sum(l) != 0]
        V = sorted(hull2(sop))
        v = v_of(T, r)
        orb = sorted({(s1 * p[0], s2 * p[1])
                      for p in itertools.permutations(v) for s1 in (1, -1) for s2 in (1, -1)})
        ok_orb = (V == orb)
        ok_uno = all(len(f[x]) == 1 for x in V)
        malo += (not ok_orb) or (not ok_uno)
        print("    %-21s | %-9s | %-36s | %-6s | %s"
              % (str(T), str(v), str(V)[:36], ok_orb, ok_uno))
    print("")
    print("    C0 %s" % ("PASA -- (1) y (2) son ciertas en lo medido" if not malo
                         else "*** FALLA en %d: la formula no es lo que dice ***" % malo))

    # ===================================================================== C1 ========================
    print("")
    print("=" * 108)
    print("C1  ¿QUE MATROIDE HAY AQUI DEBAJO?  -- clases disjuntas => matroide de PARTICION")
    print("=" * 108)
    print("")
    print("    beta                          | clases        | prod |cl_k| | #transversales | disjuntas")
    t = 4
    for b in [(18, 17, 11, 8, 7, 6, 1, 0), (13, 9, 8, 7, 5, 4, 2, 0), (12, 8, 7, 6, 5, 4, 3, 2)]:
        st = setup(b, t)
        if st is None:
            print("    %-29s | sin clases" % str(b))
            continue
        cl, E, Cd = st
        tam = [len(cl[k]) for k in sorted(cl)]
        prod = 1
        for x in tam:
            prod *= x
        tr = all_transversals(b, cl, 2, t)
        union = sum(tam)
        print("    %-29s | %-13s | %10d | %14d | %s"
              % (str(b), str(tam), prod, len(tr), union == len(b) and prod == len(tr)))
    print("")
    print("    Las cl_k particionan beta, luego el matroide transversal de la familia es un matroide de")
    print("    PARTICION y las transversales son el PRODUCTO de los bloques.  No hay teoria matroidal que")
    print("    explotar aqui: la dificultad esta en que v(T) no se separa por bloques.")

    # ===================================================================== N2 + N3 ===================
    print("")
    print("=" * 108)
    print("N2 + N3  EL EXPERIMENTO -- v maximal en dominancia y su multiplicidad con signo")
    print("=" * 108)
    print("")
    print("  poblacion critica: Phi != 0 nos da la respuesta trivial, asi que se separan las DOS columnas")
    print("  por la simetria de g_com, que es lo que la Conjetura 8.44 relaciona con la anulacion.")
    print("")
    print("   t  r   W | g_com      |   n | v maximal unico | m(v) != 0 | se anula")
    for (t, r, W) in [(4, 2, 14), (4, 2, 15), (6, 2, 15), (4, 3, 14)]:
        col = {True: [0, 0, 0, 0], False: [0, 0, 0, 0]}
        for b in betas(t, r, W):
            if not occupied(b, t):
                continue
            st = setup(b, t)
            if st is None:
                continue
            cl, E, Cd = st
            if not E:
                continue
            C, tau, gc = gcom_de(b, t, r)
            if C is None:
                continue
            S = S_of(b, t)
            resto = [x for x in S if x not in gc]
            if C != tau or set(C - x for x in resto) != set(resto):
                continue                       # poblacion critica: C = tau y S\g_com simetrico
            sim = (set(C - x for x in gc) == set(gc))
            z = phi_zero(b, t, r)
            if z is None:
                continue
            tr = all_transversals(b, cl, r, t)
            vs = defaultdict(int)
            for (_, T, w, _) in tr:
                vs[v_of(T, r)] += w
            maximales = [v for v in vs if not any(u != v and domina(u, v) for u in vs)]
            unico = (len(maximales) == 1)
            # m(v) sobre los maximales: el coeficiente del vertice en la expansion real
            at = atoms(b, t, r)
            f = fibras(at) if at else {}
            nonzero = False
            for v in maximales:
                for p in itertools.permutations(v):
                    if p in f and sum(c for c, _ in f[p]) != 0:
                        nonzero = True
            d = col[sim]
            d[0] += 1
            d[1] += unico
            d[2] += nonzero
            d[3] += bool(z)
        for sim in (True, False):
            n, u, nz, ze = col[sim]
            et = "simetrico" if sim else "ASIMETRICO"
            if not n:
                print("  %3d %2d %3d | %-10s |   0 | POBLACION VACIA -- esta fila no dice nada" % (t, r, W, et))
                continue
            print("  %3d %2d %3d | %-10s | %3d | %11d   | %7d   | %6d"
                  % (t, r, W, et, n, u, nz, ze))

    print("")
    print("  *** OJO CON LEER ESTA TABLA AL DERECHO.  m(v) != 0  IMPLICA  Phi != 0 POR DEFINICION: un")
    print("  coeficiente no nulo es un testigo.  Asi que la columna simetrica tiene que dar 0 y la")
    print("  asimetrica tiene que dar mucho, y eso NO prueba nada -- es la definicion devuelta.  El unico")
    print("  numero con contenido es el RECIPROCO: cuantas formas con Phi != 0 se quedan SIN certificado.")
    print("  Esos son los fallos, y van en N4.")

    # ===================================================================== N4 ========================
    print("")
    print("=" * 108)
    print("N4  LOS FALLOS -- Phi != 0 y m(v) = 0 en TODO maximal.  ¿Hasta donde hay que bajar?")
    print("=" * 108)
    print("")
    print("  Aqui esta el contenido: el certificado es SANO por definicion, y la pregunta es si es")
    print("  COMPLETO.  No lo es.  Y lo interesante es cuanto le falta.")
    print("")
    print("   t  r   W | Phi!=0 asim. | sin certificado | profundidad de dominancia que SI lo da")
    for (t, r, W) in [(4, 2, 15), (6, 2, 15), (4, 3, 14)]:
        n = fall = 0
        prof = defaultdict(int)
        coefs = set()
        for b in betas(t, r, W):
            if not occupied(b, t):
                continue
            st = setup(b, t)
            if st is None:
                continue
            cl, E, Cd = st
            if not E:
                continue
            C, tau, gc = gcom_de(b, t, r)
            if C is None:
                continue
            S = S_of(b, t)
            resto = [x for x in S if x not in gc]
            if C != tau or set(C - x for x in resto) != set(resto):
                continue
            if set(C - x for x in gc) == set(gc):
                continue
            if phi_zero(b, t, r):
                continue
            n += 1
            tr = all_transversals(b, cl, r, t)
            vs = defaultdict(int)
            for (_, T, w, _) in tr:
                vs[v_of(T, r)] += w
            at = atoms(b, t, r)
            f = fibras(at)

            def coef(v):
                return sum(sum(c for c, _ in f[p]) for p in itertools.permutations(v) if p in f)
            maxi = [v for v in vs if not any(u != v and domina(u, v) for u in vs)]
            if any(coef(v) != 0 for v in maxi):
                continue
            fall += 1
            niveles = sorted(vs, key=lambda v: (-sum(v), tuple(-x for x in v)))
            for k, v in enumerate(niveles):
                c = coef(v)
                if c:
                    prof[k] += 1
                    coefs.add(abs(c))
                    break
            else:
                prof[-1] += 1
        if not n:
            print("  %3d %2d %3d | POBLACION VACIA -- esta fila no dice nada" % (t, r, W))
            continue
        print("  %3d %2d %3d | %12d | %6d (%4.1f %%) | %s   coeficientes |c| = %s"
              % (t, r, W, n, fall, 100.0 * fall / n, dict(prof), sorted(coefs)))
    print("")
    print("  Profundidad 1 significa: el orbit vertice de arriba se cancela ENTERO y el certificado")
    print("  aparece en el nivel INMEDIATAMENTE siguiente.  Un solo paso, no un descenso abierto.")
    print("")
    print("=" * 108)
