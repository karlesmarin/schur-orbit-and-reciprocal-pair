# -*- coding: utf-8 -*-
# EL CRITERIO SOBRE S.  14 de agosto de 2026.
#
# layer_condition.py dio, sobre la poblacion RESTRINGIDA (ocupacion, |G|=2, C=tau), un predicado
# exacto en las cinco configuraciones: S simetrico bajo x -> tau - x.  Este guion lo ataca en la
# poblacion ENTERA, que es donde tiene que valer para ser criterio.
#
# EL ENUNCIADO CANDIDATO, escrito antes de correr:
#
#     Phi_{t,r}(lambda) == 0   <=>   falta una clase residual        [rama (a)]
#                                    o   C - S = S,  C = min S + max S
#
# donde S son los valores de las clases de EXCESO.  Notese que se dice con C y no con tau: C es
# intrinseco al beta-conjunto, tau es un dato del voraz, y si S es C-simetrico entonces min+max = C
# sale solo.  Si el criterio se sostiene, tau desaparece del enunciado.
#
# POR QUE ES LA FORMA CORRECTA, Y NO "beta simetrico".  En t=2 con e=2 se tiene S = beta y la
# condicion ES la autocomplementariedad -- el Teorema 8.6, recuperado.  Para t>=4 el conjunto S es
# mas pequeño que beta, y FRAME_RESTRICTION_KERNEL.md §6 ya habia dejado escrita la precision:
# "la parte de EXCESO de lambda es autodual", no "lambda es autodual".  Esto es esa frase, medida.
#
# LO QUE HAY QUE VIGILAR, Y POR ESO SE BARRE r=1.  En r=1 el criterio probado es la CONCENTRICIDAD,
# a1+a2 = b1+b2.  Concentrico obliga a que los dos intervalos esten ANIDADOS -- interleaved es
# imposible: a1>b1>a2>b2 con a1+a2=b1+b2 daria a1-b1 = b2-a2 < 0 -- y anidado implica S simetrico.
# Pero el reciproco NO es obvio: S = {s1>s2>s3>s4} simetrico con las clases repartidas como
# A={s1,s2}, B={s3,s4} NO es concentrico.  Si esa configuracion existe y no se anula, el candidato
# tiene un falso positivo y hay que decirlo.  Por eso r=1 entra en el barrido y no se da por hecho.
#
# CONTROLES
#   C1  r = 1 contra el Teorema 3.1, que esta PROBADO: el criterio tiene que coincidir con
#       "concentrico o falta una clase", no parecerse.
#   C2  t = 2 contra el Teorema 8.6, tambien PROBADO.
#   C3  un señuelo que TIENE que fallar: "beta entero es C-simetrico" en vez de S.  Si acertara
#       igual, este barrido no distingue las dos lecturas y no habria aprendido nada.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python criterion_S.py

import itertools
import sys

from second_stratum import setup, all_transversals
from depth_histogram import measure

CFG = [(2, 1, 16), (3, 1, 16), (4, 1, 16), (5, 1, 15), (6, 1, 15), (7, 1, 15), (8, 1, 15),
       (2, 2, 17), (2, 3, 16), (2, 4, 16), (4, 2, 18), (6, 2, 18), (8, 2, 18), (10, 2, 19),
       (12, 2, 19), (3, 2, 17), (5, 2, 17), (7, 2, 17), (9, 2, 18),
       (4, 3, 16), (6, 3, 16), (8, 3, 17), (3, 3, 16), (5, 3, 16)]


def betas(t, r, W):
    N = t + 2 * r
    for mid in itertools.combinations(range(1, W + 1), N - 1):
        yield tuple(sorted(mid, reverse=True)) + (0,)


def anat(b, t, r):
    """(Phi==0, S, C, |G|, tau) -- o None si falta una clase (rama (a), tratada aparte)."""
    st = setup(b, t)
    if st is None:
        return None
    cl, E, Cd = st
    if not E:
        return None
    S = sorted({v for k in E for v in Cd[k]})
    incr = []
    for k in E:
        c = sorted(Cd[k], reverse=True)
        incr += [c[i] + c[i + 1] for i in range(len(c) - 1)]
    incr.sort(reverse=True)
    tau = incr[r - 1] if len(incr) >= r else None
    INCR[b] = incr
    tr = all_transversals(b, cl, r, t)
    D = max(x[3] for x in tr)
    nG = sum(1 for x in tr if x[3] == D)
    zero = measure([(x[2], x[1]) for x in tr], r) is None
    return zero, S, S[0] + S[-1], nG, tau


def sym(vals, C):
    return sorted(C - v for v in vals) == sorted(vals)


def concentric(b, t):
    """el criterio PROBADO en r=1: dos clases de tamaño 2 con a1+a2 = b1+b2 (d3 = 0).
    En el perfil de tamaño tres d3 = d1+d2 > 0, luego nunca se anula."""
    cl = {}
    for v in b:
        cl.setdefault(v % t, []).append(v)
    if len(cl) < t:
        return None
    big = [sorted(v, reverse=True) for v in cl.values() if len(v) >= 2]
    if len(big) == 1:
        return False                       # perfil de tamaño tres
    A, B = big[0], big[1]
    return A[0] + A[1] == B[0] + B[1]


def run():
    print("=" * 96)
    print("CANDIDATO:  Phi == 0  <=>  falta clase  o  ( |G| = 2  y  tau = C  y  C - S = S )")
    print("=" * 96)
    print("  t  r  W    ocupadas   Phi=0   S C-sim    FP    FN   |  señuelo beta C-sim: FP  FN")
    bad = 0
    ctl_fires = 0
    global MIN, INCR
    MIN = {}
    INCR = {}
    for (t, r, W) in CFG:
        n = nz = ns = fp = fn = 0
        dfp = dfn = 0
        for b in betas(t, r, W):
            a = anat(b, t, r)
            if a is None:
                continue                    # rama (a): ambos lados dicen cero, no informa
            z, S, C, nG, tau = a
            n += 1
            nz += z
            q = sym(S, C)
            p = (nG == 2 and tau == C and q)
            ns += p
            fp += (p and not z)
            fn += (z and not p)
            inc = INCR[b]
            for nm, val in (("Q2 tau=C & S-sim", tau == C and q),
                            ("Q3 |G|=2 & S-sim", nG == 2 and q),
                            ("R1 S-sim & 2 incr = C", q and inc.count(C) >= 2),
                            ("R2 S-sim & 1 incr = C", q and inc.count(C) >= 1)):
                d = MIN.setdefault(nm, [0, 0])
                d[0] += (val and not z)
                d[1] += (z and not val)
            d = sym(b, b[0] + b[-1])        # el señuelo: beta ENTERO
            dfp += (d and not z)
            dfn += (z and not d)
        bad += fp + fn
        if dfp or dfn:
            ctl_fires += 1
        print("  %2d %2d %2d %10d %7d %8d %5d %5d   |                     %3d %3d"
              % (t, r, W, n, nz, ns, fp, fn, dfp, dfn))

    print()
    print("=" * 96)
    print("C1  r = 1 contra el Teorema 3.1 (PROBADO)")
    print("=" * 96)
    c1 = 0
    for (t, r, W) in CFG:
        if r != 1:
            continue
        agree = tot = 0
        for b in betas(t, r, W):
            a = anat(b, t, r)
            if a is None:
                continue
            z, S, C, nG, tau = a
            k = concentric(b, t)
            tot += 1
            agree += (k == z)
        c1 += tot - agree
        print("   t=%d : %d formas, el criterio probado y la expansion coinciden en %d" % (t, tot, agree))
    print("   desacuerdos: %d   <- 0 valida el instrumento contra un teorema" % c1)

    print()
    print("=" * 96)
    print("C2  t = 2 contra el Teorema 8.6 (PROBADO): S C-simetrico tiene que SER autocomplementario")
    print("=" * 96)
    c2 = 0
    for (t, r, W) in CFG:
        if t != 2:
            continue
        N = 2 + 2 * r
        same = tot = 0
        for b in betas(t, r, W):
            a = anat(b, t, r)
            if a is None:
                continue
            z, S, C, nG, tau = a
            tot += 1
            Cb = b[0] + b[-1]
            sc = (all(b[j] + b[N - 1 - j] == Cb for j in range(N)) and Cb % 2 == 0)
            same += (sc == (nG == 2 and tau == C and sym(S, C)))
        c2 += tot - same
        print("   r=%d : %d formas, las dos lecturas coinciden en %d" % (r, tot, same))
    print("   desacuerdos: %d   <- 0 dice que en t=2 el candidato ES el teorema" % c2)

    print()
    print("=" * 96)
    print("  fallos del candidato: %d      el señuelo (beta entero) falla en %d de %d configuraciones"
          % (bad, ctl_fires, len(CFG)))
    if ctl_fires == 0:
        print("  *** el señuelo no falla nunca: este barrido no separa S de beta, no mide nada")
        bad += 1
    print("  formas minimas:  " + "   ".join(
        "%s FP=%d FN=%d" % (k, v[0], v[1]) for k, v in sorted(MIN.items())))
    print("=" * 96)
    return 1 if (bad + c1 + c2) else 0


if __name__ == "__main__":
    sys.exit(run())
