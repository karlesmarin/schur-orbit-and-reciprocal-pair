# -*- coding: utf-8 -*-
# RE-MEDIDA DEL "e = t" DE LOS 12 SUPERVIVIENTES, CON M GRANDE.  13 de agosto de 2026.
#
# POR QUE.  La noche del 12, e_equals_t.py sobre las configuraciones de dim_certificate.py
# (M = 15..19) dio la contingencia
#       e < t : 216 formas objetivo,   0 supervivientes   <-- de aqui salio "e = t es NECESARIO"
#       e = t : 196 formas objetivo,  12 supervivientes
# y de ella el corolario (L2): "e = t exige t <= 2r, luego t=6 r=2, t=8 r=2, t=8 r=3 y t=10 r=2 no
# pueden dar supervivientes".  La madrugada del 13 la lectura HERMANA --"los contraejemplos de dos
# estratos tienen todos e = t"-- resulto ser un ARTEFACTO DEL RANGO: a M = 21 aparecen con e < t.
# Aqui se re-mide ESTA, que es otra medida (el segundo estrato via el certificado Delta, no el de
# abajo), con M grande y en las configuraciones donde e < t es lo unico posible.
#
# DEFINICIONES, importadas del guion original y no reescritas: setup / all_transversals / inv_of de
# second_stratum.py, dim_gl / halves de dim_certificate.py, P_poly de second_vanishes.py.
#
#   superviviente  :=  forma de la poblacion objetivo ([Phi]_top = 0 y (i) falsa) con Delta = 0,
#                      Delta = sum_{g in G2} w(g) dim(atil) dim(astar).
#
# COLUMNAS
#   C0  ACEPTACION, fatal: con las configuraciones ORIGINALES hay que reproducir 412 objetivo,
#       12 supervivientes y la contingencia 216/0 + 196/12.  Si no salen, mido otra cosa.
#   N1  la misma contingencia con M grande, por configuracion.
#   N2  (L1) supervivientes con e < t.  Si aparece uno, e = t deja de ser necesario.
#   N3  (L2) las cuatro configuraciones declaradas imposibles: t=6 r=2, t=8 r=2, t=8 r=3, t=10 r=2.
#   N4  Delta = 0  vs  [Phi]_{D2} = 0 de verdad (el polinomio).  Delta != 0 => no nulo es SANO por
#       construccion; que Delta = 0 implique nulo es una MEDIDA (132/132 en el rango viejo).  Un
#       superviviente con el polinomio no nulo NO es un superviviente.
#   N5  cada superviviente, con la expansion ENTERA: Phi_t != 0 y en que estrato aparece.  La noche
#       del 12 dijo "siempre en D1 - 4".
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python survivors_wide.py

import itertools
import sys
from collections import defaultdict, Counter

from second_stratum import setup, all_transversals, inv_of, perm_sign
from second_vanishes import P_poly
from dim_certificate import dim_gl, halves, CONFIGS as CONF_OLD

# donde e < t es lo UNICO posible (t > 2r) van los M mas grandes que aguanta el reloj
CONF_NEW = [(6, 2, 19), (8, 2, 20), (8, 3, 20), (8, 3, 21), (10, 2, 21),
            (4, 2, 19), (6, 3, 19), (6, 4, 20)]


def full_expansion(tm, r):
    n = 2 * r
    D = defaultdict(int)
    for w, T in tm:
        for q in itertools.permutations(range(n)):
            e = [0] * r
            for a in range(n):
                e[q[a] // 2] += (T[a] if q[a] % 2 == 0 else -T[a])
            D[tuple(e)] += w * perm_sign(list(q))
    return {k: v for k, v in D.items() if v}


def scan(t, r, M, keep=True):
    """(n_objetivo, contingencia, supervivientes) con la definicion del guion original."""
    N = t + 2 * r
    if M < N - 1:
        return None
    cont = Counter()                     # (e==t, superviviente) -> formas
    surv = []
    for comb in itertools.combinations(range(M + 1), N):
        beta = tuple(sorted(comb, reverse=True))
        st = setup(beta, t)
        if st is None:
            continue
        cl, E, Cd = st
        if not E:
            continue
        tr = all_transversals(beta, cl, r, t)
        D = max(x[3] for x in tr)
        G = [x for x in tr if x[3] == D]
        if len(G) != 2:
            continue
        a, b = G
        if not (inv_of(a[1], r) == inv_of(b[1], r) and a[2] == -b[2]):
            continue                                  # [Phi]_top != 0
        S = sorted({v for k in E for v in Cd[k]})
        C = S[0] + S[-1]
        if set(C - v for v in S) == set(S):
            continue                                  # (i) cierta
        rest = [x for x in tr if x[3] < D]
        if not rest:
            continue
        D2 = max(x[3] for x in rest)
        G2 = [x for x in rest if x[3] == D2]
        delta = 0
        for (_, T, w, _) in G2:
            at, ast = halves(T, r)
            delta += w * dim_gl(at) * dim_gl(ast)
        is_surv = (delta == 0)
        cont[(len(E) == t, is_surv)] += 1
        if is_surv and keep:
            acc = defaultdict(int)
            for (_, T, w, _) in G2:
                for k, v in P_poly(T, r).items():
                    acc[k] += w * v
            surv.append(dict(t=t, r=r, beta=beta, e=len(E), nG2=len(G2), D1=D,
                             poly_zero=not any(acc.values()),
                             tm=[(x[2], x[1]) for x in tr]))
    return sum(cont.values()), cont, surv


def contingency(cont):
    return (cont[(False, False)], cont[(False, True)],
            cont[(True, False)], cont[(True, True)])


if __name__ == "__main__":
    # ================================================================= C0 ============================
    print("=" * 104)
    print("C0  ACEPTACION: reproducir e_equals_t.py en SUS configuraciones (objetivo 412, surv 12,")
    print("    contingencia  e<t: 216 formas / 0 supervivientes   e=t: 196 formas / 12)")
    print("=" * 104)
    print("")
    print("     t   r    M | objetivo | e<t surv | e=t surv")
    print("  " + "-" * 88)
    tot = 0
    CC = Counter()
    old_surv = []
    for (t, r, M) in CONF_OLD:
        out = scan(t, r, M)
        if out is None:
            continue
        n, cont, sv = out
        lf, lt, ef, et = contingency(cont)
        print("  %4d %3d %4d | %8d | %3d %4d | %3d %4d" % (t, r, M, n, lf + lt, lt, ef + et, et))
        tot += n
        CC.update(cont)
        old_surv += sv
        sys.stdout.flush()
    lf, lt, ef, et = contingency(CC)
    print("")
    print("  objetivo %d (esperado 412) | e<t: %d formas, %d supervivientes (esperado 216, 0)"
          % (tot, lf + lt, lt))
    print("  supervivientes %d (esperado 12) | e=t: %d formas, %d supervivientes (esperado 196, 12)"
          % (lt + et, ef + et, et))
    ok = (tot == 412 and lt + et == 12 and lf + lt == 216 and lt == 0 and ef + et == 196 and et == 12)
    print("")
    print("  C0 %s" % ("PASA -- son sus numeros exactos." if ok else
                       "FALLA -- estoy midiendo otra cosa, el resto NO vale."))
    if not ok:
        print("DONE (veredicto suspendido)")
        raise SystemExit(1)

    # ================================================================= N1-N3 =========================
    print("")
    print("=" * 104)
    print("N1-N3  la misma contingencia con M GRANDE")
    print("=" * 104)
    print("")
    print("     t   r    M    N  t>2r | objetivo | e<t: formas  surv | e=t: formas  surv")
    print("  " + "-" * 100)
    NEW = Counter()
    new_surv = []
    for (t, r, M) in CONF_NEW:
        out = scan(t, r, M)
        if out is None:
            print("  %4d %3d %4d %4d | SALTADA (M < N-1)" % (t, r, M, t + 2 * r))
            continue
        n, cont, sv = out
        lf, lt, ef, et = contingency(cont)
        print("  %4d %3d %4d %4d %5s | %8d | %11d %5d | %11d %5d"
              % (t, r, M, t + 2 * r, "si" if t > 2 * r else "no", n, lf + lt, lt, ef + et, et))
        NEW.update(cont)
        new_surv += sv
        sys.stdout.flush()

    lf, lt, ef, et = contingency(NEW)
    print("")
    print("  N2 (L1) 'e = t es NECESARIO entre los supervivientes':")
    print("       e < t : %d formas objetivo, %d supervivientes" % (lf + lt, lt))
    print("       e = t : %d formas objetivo, %d supervivientes" % (ef + et, et))
    print("       %s" % ("(L1) SOBREVIVE en este rango." if lt == 0 else
                         "*** (L1) REFUTADA: hay supervivientes con e < t ***"))
    print("")
    print("  N3 (L2) las configuraciones declaradas imposibles por 'e = t exige t <= 2r':")
    for key in [(6, 2), (8, 2), (8, 3), (10, 2)]:
        k = [s for s in new_surv if (s['t'], s['r']) == key]
        print("       t=%d r=%d : %d supervivientes%s"
              % (key[0], key[1], len(k), "" if not k else "   *** (L2) REFUTADA ***"))

    # ================================================================= N4, N5 ========================
    print("")
    print("=" * 104)
    print("N4/N5  cada superviviente: Delta = 0 contra el POLINOMIO, y donde aparece Phi_t")
    print("=" * 104)
    print("")
    print("     t   r  e   t  |G2|  beta                                     poly=0   D1   1er  prof.")
    print("  " + "-" * 100)
    prof = Counter()
    n_poly_bad = n_zero = 0
    for s in sorted(new_surv, key=lambda x: (x['t'], x['r'], x['beta']))[:40]:
        FE = full_expansion(s['tm'], s['r'])
        if not s['poly_zero']:
            n_poly_bad += 1
        if not FE:
            n_zero += 1
            line = "  Phi_t == 0 (!!)"
            mx = d = 0
        else:
            mx = max(sum(k) for k in FE)
            d = s['D1'] - mx
            prof[d] += 1
            line = "%5d %5d" % (mx, d)
        print("  %4d %3d %2d %3d %5d  %-40s %6s %4d %s"
              % (s['t'], s['r'], s['e'], s['t'], s['nG2'], str(list(s['beta'])),
                 "si" if s['poly_zero'] else "NO", s['D1'], line))
        sys.stdout.flush()
    if len(new_surv) > 40:
        print("  ... y %d mas" % (len(new_surv) - 40))
    print("")
    print("  N4  supervivientes con Delta = 0 pero POLINOMIO no nulo (no son supervivientes): %d de %d"
          % (n_poly_bad, len(new_surv)))
    print("  N5  profundidad del primer estrato no nulo (la noche del 12: SIEMPRE D1 - 4): %s"
          % dict(prof))
    print("      supervivientes con Phi_t == 0 (contraejemplo a la necesidad de (i)): %d" % n_zero)
    print("")
    print("DONE")
