# -*- coding: utf-8 -*-
# EL TEOREMA 2 DE PRASAD (Israel J. Math 211, 2016), VERIFICADO --- Y PUESTO AL LADO DEL NUESTRO.
# 25 de agosto de 2026.
#
# POR QUE.  Prasad pregunta si hay teoria para las POTENCIAS de a_P.  Resulta que el escribio la
# teoria para a_Q (Coxeter) en 2016, y su enunciado tiene EXACTAMENTE la forma de lo que hemos
# medido hoy en tipo C.  Antes de decirselo hay que (1) verificar SU teorema con nuestra maquina,
# y (2) poner los dos enunciados al lado y comprobar que el diccionario es el que decimos.
#
# SU TEOREMA 2, literal del PDF (_papers/_prasad_israe.txt):
#   "Theta(t.c_n) is not identically zero if and only if the mn-tuple of integers appearing in
#    lambda + rho_mn represents each residue class in Z/n exactly m-times.  Assume this to be the
#    case.  [...]  Theta_pi(t.c_n) = +- Theta_1(t^n) Theta_2(t^n) ... Theta_n(t^n)."
# con mu_i + rho_m = [lambda + rho_mn - i]/n,  y Theta_i el caracter de GL_m.
#
# EN t = 1 el alfabeto es mu_n con multiplicidad m, luego
#     H(z) = prod_{w^n=1} 1/(1-wz)^m = 1/(1-z^n)^m   ->   h_k = C(j+m-1, m-1) si k = jn, 0 si no.
# Enteros: determinante exacto.  Y Theta_i(1) = dim, por la formula de dimensiones de Weyl de GL_m.
#
# QUE SE MIDE
#   C0  CONTROL: la Jacobi-Trudi para GL contra el bialternante, en un punto GENERICO.
#   P1  SU criterio de anulacion: chi != 0  <=>  cada residuo mod n exactamente m veces.
#   P2  SU formula del valor en t=1: |chi| = producto de las dimensiones de los n bloques GL_m.
#       CONTROL: un senuelo que use m-1 o m+1 bloques debe FALLAR.
#   D   el DICCIONARIO con lo nuestro, impreso para leerlo de un vistazo.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python prasad2016_vs_nuestro.py > prasad2016_vs_nuestro_OUT.txt 2>&1

import itertools
import random
import sys
from collections import Counter
from math import comb

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SEP = "=" * 94


def det_int(M):
    n = len(M)
    if n == 0:
        return 1
    M = [row[:] for row in M]
    sign, prev = 1, 1
    for c in range(n - 1):
        if M[c][c] == 0:
            piv = next((r for r in range(c + 1, n) if M[r][c] != 0), None)
            if piv is None:
                return 0
            M[c], M[piv] = M[piv], M[c]
            sign = -sign
        for r in range(c + 1, n):
            for k in range(c + 1, n):
                M[r][k] = (M[r][k] * M[c][c] - M[r][c] * M[c][k]) // prev
            M[r][c] = 0
        prev = M[c][c]
    return sign * M[n - 1][n - 1]


def det_f(M):
    n = len(M)
    M = [row[:] for row in M]
    d = 1.0
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        if abs(M[p][c]) < 1e-12:
            return 0.0
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            for k in range(c, n):
                M[r][k] -= f * M[c][k]
    return d


def schur_JT(lam, H, L):
    """s_lambda = det(h_{lam_i - i + j})_{1<=i,j<=L}."""
    return [[H(lam[i] - (i + 1) + (j + 1)) for j in range(L)] for i in range(L)]


def dim_GL(mu):
    """Formula de dimensiones de Weyl para GL_m:  prod_{k<l} (b_k-b_l)/(l-k),  b = mu + rho_m."""
    m = len(mu)
    b = [mu[i] + (m - 1 - i) for i in range(m)]
    num = den = 1
    for k in range(m):
        for l in range(k + 1, m):
            num *= (b[k] - b[l])
            den *= (l - k)
    return num // den


# ------------------------------------------------------------------ C0
print(SEP)
print("C0 -- CONTROL: la Jacobi-Trudi de GL reproduce el bialternante en un punto GENERICO")
random.seed(11)
for N in (3, 4):
    xs = [random.uniform(1.2, 2.4) for _ in range(N)]
    K = 40
    h = [0.0] * (K + 1)
    h[0] = 1.0
    for x in xs:
        acc, new = 0.0, [0.0] * (K + 1)
        for k in range(K + 1):
            acc = acc * x + h[k]
            new[k] = acc
        h = new
    Hf = lambda k: (h[k] if 0 <= k < len(h) else 0.0)
    ok = tot = 0
    for lam in itertools.product(range(4), repeat=N):
        if any(lam[i] < lam[i + 1] for i in range(N - 1)):
            continue
        tot += 1
        e = [lam[i] + (N - 1 - i) for i in range(N)]
        num = det_f([[x ** ei for x in xs] for ei in e])
        den = det_f([[x ** (N - 1 - i) for x in xs] for i in range(N)])
        ref = num / den
        val = det_f(schur_JT(list(lam), Hf, N))
        ok += abs(val - ref) <= 1e-6 * max(1.0, abs(ref))
    print(f"     N={N}: {ok}/{tot}")

# ------------------------------------------------------------------ P1, P2
print(SEP)
print("P1/P2 -- EL TEOREMA 2 DE PRASAD (2016), en t=1.   alfabeto = mu_n con multiplicidad m")
print("         h_k = C(k/n + m-1, m-1) si n | k, y 0 si no.")
print()
print(f"{'m':>3} {'n':>3} {'mn':>4} {'pesos':>7} {'P1 criterio':>13} {'P2 valor':>11} "
      f"{'senuelo n-1 bloques':>20} {'senuelo n+1':>13}")
for (m, n, top) in [(1, 3, 6), (2, 2, 6), (2, 3, 5), (3, 2, 5), (2, 4, 4), (3, 3, 3), (4, 2, 4)]:
    N = m * n
    Hc = lambda k: (comb(k // n + m - 1, m - 1) if (k >= 0 and k % n == 0) else 0)
    ok1 = ok2 = tot = d1 = d2 = 0
    for lam in itertools.product(range(top + 1), repeat=N):
        if any(lam[i] < lam[i + 1] for i in range(N - 1)):
            continue
        tot += 1
        val = det_int(schur_JT(list(lam), Hc, N))
        a = [lam[i] + (N - 1 - i) for i in range(N)]          # lambda + rho_mn
        cnt = Counter(x % n for x in a)
        crit = all(cnt.get(i, 0) == m for i in range(n))
        ok1 += (crit == (val != 0))
        if crit:
            # los n bloques:  mu_i + rho_m = [a - i]/n  sobre los a con residuo i
            pred = 1
            for i in range(n):
                b = sorted((x - i) // n for x in a if x % n == i)[::-1]
                mu = [b[k] - (m - 1 - k) for k in range(m)]
                pred *= dim_GL(mu)
            ok2 += (abs(val) == pred)
            # senuelos: quitar un bloque / duplicar uno
            bl = []
            for i in range(n):
                b = sorted((x - i) // n for x in a if x % n == i)[::-1]
                bl.append(dim_GL([b[k] - (m - 1 - k) for k in range(m)]))
            p_menos = 1
            for v in bl[:-1]:
                p_menos *= v
            p_mas = pred * (bl[0] if bl else 1)
            d1 += (abs(val) == p_menos)
            d2 += (abs(val) == p_mas)
        else:
            ok2 += 1
    nz = sum(1 for _ in [0])
    print(f"{m:>3} {n:>3} {N:>4} {tot:>7} {str(ok1)+'/'+str(tot):>13} {str(ok2)+'/'+str(tot):>11} "
          f"{d1:>20} {d2:>13}")
print()
print("     (Las dos ultimas columnas son SENUELOS: cuantas veces acertaria una formula con un")
print("      bloque de menos o uno de mas.  Cuanto mas bajas, mas discrimina el test.)")

# ------------------------------------------------------------------ D
print(SEP)
print("D -- EL DICCIONARIO.   Su Teorema 2 (tipo A, a_Q)  <->  lo nuestro (tipo C, a_P)")
print()
filas = [
    ("el elemento", "t . c_n,  c_n = clase de Coxeter", "a_P^k  (segundo elemento de Kostant)"),
    ("el orden", "n", "q = t/d,  t = 2m+2,  d = gcd(k,t)"),
    ("la multiplicidad", "m", "d"),
    ("las clases", "residuos mod n, SIN plegar", "clases PLEGADAS mod q:  cls(x)=min(x, -x)"),
    ("criterio de anulacion", "cada residuo EXACTAMENTE m veces", "cada clase A LO SUMO d veces"),
    ("", "(no hay clases fijas)", "y a lo sumo floor(d/2) en las FIJAS 0 y q/2"),
    ("el valor", "+- prod de n caracteres de GL_m", "+- prod de dimensiones de Weyl de R_P(q)"),
    ("en t=1 / en el elemento", "+- prod de n DIMENSIONES", "+- prod de dimensiones"),
    ("el caso degenerado", "m = 1  ->  {0,+-1}  = KOSTANT", "d = 1  ->  {0,+-1}  = a_P, Kostant"),
]
w = max(len(f[0]) for f in filas)
print(f"  {'':<{w}} | {'PRASAD 2016, tipo A':<40} | lo nuestro, tipo C")
print(f"  {'-'*w}-+-{'-'*40}-+-{'-'*44}")
for a, b, c in filas:
    print(f"  {a:<{w}} | {b:<40} | {c}")
print()
print("  LA DIFERENCIA, en una frase:  el PLEGADO.  En tipo A los residuos no se pliegan y la")
print("  condicion es una IGUALDAD (exactamente m); en tipo C se pliegan por r -> q-r, la condicion")
print("  se afloja a una DESIGUALDAD (a lo sumo d), y las dos clases fijas 0 y q/2 --- las que la")
print("  involucion no mueve --- llevan la mitad de capacidad.  Esa asimetria ES la raiz larga.")

# ------------------------------------------------------------------ W
print(SEP)
print("W -- POR QUE una es IGUALDAD y la otra DESIGUALDAD.  No se postula: se cuenta.")
print()
print("   TIPO A:  hay mn beta-numeros y n clases de capacidad m.  Capacidad total = mn = numero de")
print("            filas.  SATURA.  Luego 'a lo sumo m' y 'exactamente m' son lo mismo, y el")
print("            enunciado natural es la IGUALDAD.")
print()
print("   TIPO C:  hay m = dq/2 - 1 filas.  Las clases plegadas mod q son floor(q/2)+1, de las que")
print("            son FIJAS una (q impar) o dos (q par).  Capacidad total = d*(no fijas) +")
print("            floor(d/2)*(fijas), y sale")
print("                  d PAR   ->  capacidad = dq/2      ->  holgura = 1")
print("                  d IMPAR ->  capacidad = dq/2 - 1  ->  holgura = 0  (satura, como en A)")
print()
from math import gcd as _gcd
mal = tot = 0
ejemplos = []
for m_ in range(2, 60):
    t_ = 2 * m_ + 2
    for k_ in range(1, 40):
        d_ = _gcd(k_, t_)
        q_ = t_ // d_
        if q_ < 3:
            continue
        fij = 1 + (1 if q_ % 2 == 0 else 0)
        nof = (q_ // 2 - 1) if q_ % 2 == 0 else (q_ - 1) // 2
        cap = d_ * nof + (d_ // 2) * fij
        tot += 1
        if cap - m_ != (1 if d_ % 2 == 0 else 0):
            mal += 1
        if m_ <= 7 and k_ <= 4:
            ejemplos.append((m_, k_, d_, q_, m_, cap, cap - m_))
print(f"   {'m':>3} {'k':>3} {'d':>3} {'q':>4} {'filas':>6} {'capacidad':>10} {'holgura':>8}")
for e in ejemplos[:16]:
    print(f"   {e[0]:>3} {e[1]:>3} {e[2]:>3} {e[3]:>4} {e[4]:>6} {e[5]:>10} {e[6]:>8}")
print(f"\n   Comprobado en {tot} casos (m=2..59, k=1..39): {mal} fallos.")
print("   CONSECUENCIA:  cuando d es IMPAR el criterio de tipo C tambien satura, y se puede leer")
print("   como igualdad --- o sea que es literalmente el enunciado de Prasad 2016 plegado.  La")
print("   desigualdad solo hace falta cuando d es PAR, y la unidad de holgura viene ENTERA de las")
print("   clases fijas: cada una lleva floor(d/2) en vez de d, y d - 2*floor(d/2) = d mod 2.")
print(SEP)
