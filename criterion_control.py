# -*- coding: utf-8 -*-
"""CONTROL INDEPENDIENTE del criterio del 14 de agosto, en Python exacto y sin Sage.

POR QUE EXISTE.  criterion_S.py sale entero de un instrumento: setup(), all_transversals() y
measure() -- Laplace, estratos y voraz.  Un fallo ahi convertiria 24 configuraciones limpias en 24
acuerdos consigo mismo.  Este fichero no comparte una linea con aquel.

LA VIA.  Jacobi-Trudi.  El alfabeto es mu_t union {z_j, 1/z_j}, y su serie generatriz de las h
colapsa, porque la orbita entera se contrae:

        prod_{k<t} (1 - zeta^k u) = 1 - u^t                 (las raices t-esimas de la unidad)

    =>  H(u) = 1/(1 - u^t) * prod_j 1/((1 - z_j u)(1 - u/z_j)),      s_lambda = det( h_{L_i - i + j} ).

No aparece ninguna raiz de la unidad: el calculo vive sobre Q.  Ni bialternante, ni expansion de
Laplace, ni estratificacion por grados, ni greedy.

QUE PRUEBA ESTE CONTROL, EXACTAMENTE -- y la asimetria es lo importante.  Se evalua en puntos
RACIONALES, con Fraction, sin coma flotante:

  * si algun punto da un valor != 0, entonces Phi NO se anula identicamente.  Eso es una PRUEBA.
    Luego este control REFUTA EXACTAMENTE los falsos positivos del criterio: si el criterio dice
    "se anula" y no se anula, un solo punto lo caza.
  * si todos los puntos dan 0, es evidencia fuerte y no prueba.  Se usan 12 puntos independientes
    y se dice asi.

O sea: la direccion que el criterio podria estar inventando -- ceros que no existen -- queda cerrada
sin margen; la otra queda medida.

CONTROLES DE ESTE FICHERO
  K1  en r = 1 el criterio PROBADO es la concentricidad (Teorema 3.1).  El instrumento tiene que
      reproducirlo, o no es un instrumento.
  K2  un señuelo que TIENE que fallar: pedir C - beta = beta (beta entero) en vez de C - S = S.
  K3  el cedazo tiene que ser capaz de dar positivos: se cuenta en cuantas formas encuentra un
      punto no nulo.  Si fuera 0, no estaria midiendo.

CHECKPOINT.  Cada forma se escribe a criterion_control_CKPT.jsonl con flush; al arrancar se lee y
se salta lo hecho.  Progreso por configuracion con flush.

Authors: Carles Marin, Claude (AI assistant).
"""
import itertools
import json
import os
import random
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
CKPT = os.path.join(HERE, "criterion_control_CKPT.jsonl")
NPTS = 12
CFG = [(2, 1, 14), (3, 1, 14), (4, 1, 14), (5, 1, 13), (6, 1, 13),
       (2, 2, 15), (4, 2, 16), (6, 2, 16), (8, 2, 17), (3, 2, 15), (5, 2, 15),
       (2, 3, 14), (4, 3, 14), (6, 3, 15)]


def betas(t, r, W):
    N = t + 2 * r
    for mid in itertools.combinations(range(1, W + 1), N - 1):
        yield tuple(sorted(mid, reverse=True)) + (0,)


def lam_of(b):
    N = len(b)
    return [b[i] - (N - 1 - i) for i in range(N)]


def hseries(t, r, zs, D):
    """h_0..h_D del alfabeto, por convolucion de las series exactas.  Todo en Fraction."""
    h = [F(0)] * (D + 1)
    for k in range(0, D + 1, t):          # 1/(1-u^t)
        h[k] = F(1)
    for z in zs:
        for a in (z, 1 / z):              # 1/(1-a u) = sum a^k u^k
            acc = [F(0)] * (D + 1)
            run = F(0)
            pw = F(1)
            # convolucion con la serie geometrica, por recurrencia: acc_k = h_k + a*acc_{k-1}
            for k in range(D + 1):
                run = h[k] + a * run
                acc[k] = run
            h = acc
    return h


def det_frac(M):
    n = len(M)
    M = [row[:] for row in M]
    d = F(1)
    for c in range(n):
        piv = next((rr for rr in range(c, n) if M[rr][c] != 0), None)
        if piv is None:
            return F(0)
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            d = -d
        d *= M[c][c]
        pv = M[c][c]
        for rr in range(c + 1, n):
            if M[rr][c] == 0:
                continue
            f = M[rr][c] / pv
            for k in range(c, n):
                M[rr][k] -= f * M[c][k]
    return d


def value(lam, t, r, zs):
    """s_lambda en el alfabeto, exacto, por Jacobi-Trudi."""
    N = t + 2 * r
    L = list(lam) + [0] * (N - len(lam))
    D = L[0] + N + 1
    h = hseries(t, r, zs, D)
    def H(k):
        return F(0) if k < 0 or k > D else h[k]
    return det_frac([[H(L[i] - i + j) for j in range(N)] for i in range(N)])


def excess(b, t):
    cl = {}
    for v in b:
        cl.setdefault(v % t, []).append(v)
    if len(cl) < t:
        return None
    S, incr = [], []
    for k, vs in cl.items():
        if len(vs) >= 2:
            c = sorted(vs, reverse=True)
            S += c
            incr += [c[i] + c[i + 1] for i in range(len(c) - 1)]
    if not S:
        return None
    S = sorted(S)
    return S, S[0] + S[-1], incr


def concentric(b, t):
    cl = {}
    for v in b:
        cl.setdefault(v % t, []).append(v)
    big = [sorted(v, reverse=True) for v in cl.values() if len(v) >= 2]
    if len(big) == 1:
        return False
    return big[0][0] + big[0][1] == big[1][0] + big[1][1]


def main():
    random.seed(20260814)
    done = {}
    if os.path.exists(CKPT):
        with open(CKPT) as fh:
            for line in fh:
                line = line.strip()
                if line:
                    d = json.loads(line)
                    done[d["k"]] = d["z"]
        print("checkpoint: %d formas ya calculadas" % len(done), flush=True)
    out = open(CKPT, "a")

    print("=" * 100, flush=True)
    print("CONTROL INDEPENDIENTE -- Jacobi-Trudi exacto sobre Q.  Sin Laplace, sin estratos, sin voraz.", flush=True)
    print("=" * 100, flush=True)
    print("  t  r  W    formas   nulos  criterio    FP*   FN  | K1 r=1 | K2 señuelo FP FN | K3 no-nulos",
          flush=True)
    bad = k1bad = decoy = 0
    for (t, r, W) in CFG:
        n = nz = ncrit = fp = fn = dfp = dfn = k1 = nonz = 0
        for b in betas(t, r, W):
            e = excess(b, t)
            if e is None:
                continue
            S, C, incr = e
            key = "%d,%d,%s" % (t, r, ",".join(map(str, b)))
            if key in done:
                z = bool(done[key])
            else:
                lam = lam_of(b)
                z = True
                for _ in range(NPTS):
                    pts = [F(random.randint(3, 400), random.randint(3, 400)) for _ in range(r)]
                    if len(set(pts)) < r:
                        continue
                    if value(lam, t, r, pts) != 0:
                        z = False
                        break
                out.write(json.dumps({"k": key, "z": int(z)}) + "\n")
                out.flush()
            n += 1
            nz += z
            nonz += (not z)
            crit = sorted(C - v for v in S) == S and (C in incr)
            ncrit += crit
            fp += (crit and not z)         # <- REFUTACION EXACTA
            fn += (z and not crit)
            d = sorted(C - v for v in b) == sorted(b)
            dfp += (d and not z)
            dfn += (z and not d)
            if r == 1:
                k1 += (concentric(b, t) != z)
        bad += fp + fn
        k1bad += k1
        if dfp or dfn:
            decoy += 1
        print("  %2d %2d %2d %9d %7d %9d %6d %4d  | %6s | %14d %2d | %11d"
              % (t, r, W, n, nz, ncrit, fp, fn, ("%d" % k1) if r == 1 else "-", dfp, dfn, nonz),
              flush=True)
    out.close()
    print(flush=True)
    print("  FP* son refutaciones EXACTAS: un punto racional no nulo prueba que no se anula.", flush=True)
    print("  fallos del criterio: %d      K1 (contra el Teorema 3.1): %d desacuerdos" % (bad, k1bad),
          flush=True)
    print("  K2 el señuelo falla en %d de %d configuraciones" % (decoy, len(CFG)), flush=True)
    if decoy == 0:
        print("  *** el señuelo no falla: el control no separa S de beta", flush=True)
        bad += 1
    print("VEREDICTO: %s" % ("CONFIRMA" if bad + k1bad == 0 else "DISCREPA"), flush=True)
    return 1 if bad + k1bad else 0


if __name__ == "__main__":
    sys.exit(main())
