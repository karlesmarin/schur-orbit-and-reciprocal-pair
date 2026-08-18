# -*- coding: utf-8 -*-
# POR QUE Delta != 0.  De medir el fenomeno a buscarle el mecanismo.
#
# Delta_D = sum_i (+-1) * dim(atil_i) * dim(astar_i):  signos +-1 y ENTEROS POSITIVOS.  Una suma asi
# solo se anula si los dos lados se equilibran EXACTAMENTE.  Tres mecanismos lo impedirian, y los
# tres son PRUEBA -- no medida -- cuando se dan:
#
#   H-SIG   SIGNO.  Si todos los terminos tienen el MISMO signo, Delta != 0.  El mas elemental de
#           todos, y se me paso en la primera pasada: el residuo de 2 terminos lo delato.  Con dos
#           terminos de signo opuesto, Delta != 0 equivale a que las dimensiones difieran, o sea a
#           la dominancia; si la dominancia falla con dos terminos es porque comparten signo.
#   H-DOM   DOMINANCIA.  Si  max_i |d_i|  >  suma de los otros |d_j|,  entonces Delta != 0.
#   H-VAL   VALUACION 2-ADICA.  Si un solo i alcanza  min_i v_2(d_i),  entonces Delta != 0,
#           porque al dividir por 2^min queda exactamente un sumando impar.
#
# Ninguno puede valer siempre: en las formas donde Delta = 0 los TRES tienen que fallar.  Esa es la
# tabla de contingencia, y es el control -- si alguno "acertara" tambien alli, estaria mal calculado.
#
# LO QUE SE MIDE, sobre los TRES primeros estratos de cada forma
#   K0  contingencia (Delta != 0) x (todos el mismo signo).
#   K1  contingencia (Delta != 0) x (dominancia).   K2  contingencia (Delta != 0) x (v2 unica).
#   K3  cobertura conjunta: en que fraccion de los estratos con Delta != 0 dispara ALGUNO de los tres.
#       Eso es lo que se podria convertir en teorema; el resto queda como residuo explicito.
#   K4  el residuo: estratos con Delta != 0 donde NINGUNO dispara.  Se listan sus |G_k| y su numero
#       de terminos, que es por donde habria que seguir.
#   K5  no vacuidad: n de estratos con un solo termino (ahi la dominancia es trivial y NO cuenta
#       como mecanismo; se separan).
#   K6  SENUELO: valuacion 3-adica en vez de 2-adica.  No hay razon para que el 3 sea especial, asi
#       que deberia cubrir MENOS.  Si cubriera lo mismo, es que la valuacion no esta capturando
#       nada estructural y H-VAL seria una coincidencia de rango.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python why_delta.py

import itertools
from collections import Counter

from second_stratum import setup, all_transversals
from depth import dims

CONFIGS = [(4, 2, 15), (4, 3, 15), (6, 2, 16), (6, 3, 17), (8, 2, 17), (6, 4, 18)]


def v_p(n, p):
    if n == 0:
        return 10 ** 9
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def strata_terms(beta, t, r, nstrata=3):
    """[(grado, [terminos enteros con signo])] para los nstrata grados mas altos."""
    st = setup(beta, t)
    if st is None:
        return None
    cl, E, Cd = st
    if not E:
        return None
    tr = all_transversals(beta, cl, r, t)
    idx = range(2 * r)
    by = {}
    for (_, T, w, _) in tr:
        sT = sum(T)
        for R in itertools.combinations(idx, r):
            A = tuple(T[a] for a in R)
            B = tuple(T[a] for a in idx if a not in R)
            D = 2 * sum(A) - sT
            eps = -1 if (sum(R) % 2) else 1
            by.setdefault(D, []).append(w * eps * dims(A, B, r))
    return [(D, by[D]) for D in sorted(by, reverse=True)[:nstrata]]


def main():
    k0 = Counter()
    k1 = Counter()
    k2 = Counter()
    k6 = Counter()
    single = 0
    joint_hit = joint_n = 0
    residue = Counter()
    tot = 0

    for (t, r, M) in CONFIGS:
        N = t + 2 * r
        for comb in itertools.combinations(range(M + 1), N):
            beta = tuple(sorted(comb, reverse=True))
            S = strata_terms(beta, t, r)
            if S is None:
                continue
            tot += 1
            for (D, terms) in S:
                d = sum(terms)
                if len(terms) == 1:
                    single += 1
                    continue                    # dominancia trivial: no cuenta como mecanismo
                a = sorted((abs(x) for x in terms), reverse=True)
                dom = a[0] > sum(a[1:])
                sig = all(x > 0 for x in terms) or all(x < 0 for x in terms)
                vs = [v_p(abs(x), 2) for x in terms]
                mn = min(vs)
                val = vs.count(mn) == 1
                vs3 = [v_p(abs(x), 3) for x in terms]
                mn3 = min(vs3)
                val3 = vs3.count(mn3) == 1
                k0[(d != 0, sig)] += 1
                k1[(d != 0, dom)] += 1
                k2[(d != 0, val)] += 1
                k6[(d != 0, val3)] += 1
                if d != 0:
                    joint_n += 1
                    if sig or dom or val:
                        joint_hit += 1
                    else:
                        residue[len(terms)] += 1
        print("   hecho t=%d r=%d M=%d (acumulado %d formas)" % (t, r, M, tot), flush=True)

    print("")
    print("K5 estratos con UN SOLO termino (dominancia trivial, excluidos) : %d" % single)
    print("")
    print("K0 contingencia  (Delta != 0)  x  TODOS EL MISMO SIGNO:")
    for k in sorted(k0):
        print("      Delta!=0=%-5s mismosigno=%-5s : %d" % (k[0], k[1], k0[k]))
    print("")
    print("K1 contingencia  (Delta != 0)  x  DOMINANCIA:")
    for k in sorted(k1):
        print("      Delta!=0=%-5s dominancia=%-5s : %d" % (k[0], k[1], k1[k]))
    print("")
    print("K2 contingencia  (Delta != 0)  x  v_2 MINIMA UNICA:")
    for k in sorted(k2):
        print("      Delta!=0=%-5s v2unica=%-5s : %d" % (k[0], k[1], k2[k]))
    print("")
    print("K6 SENUELO  v_3 minima unica (deberia cubrir MENOS):")
    for k in sorted(k6):
        print("      Delta!=0=%-5s v3unica=%-5s : %d" % (k[0], k[1], k6[k]))
    print("")
    if joint_n:
        print("K3 COBERTURA CONJUNTA: alguno de los TRES dispara en %d de %d estratos con Delta != 0"
              " (%.2f %%)" % (joint_hit, joint_n, 100.0 * joint_hit / joint_n))
    print("K4 RESIDUO (Delta != 0 y NINGUNO de los tres dispara), por numero de terminos:")
    for k in sorted(residue):
        print("      %3d terminos : %d" % (k, residue[k]))
    print("      total residuo : %d" % sum(residue.values()))


if __name__ == "__main__":
    main()
