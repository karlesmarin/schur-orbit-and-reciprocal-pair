# -*- coding: utf-8 -*-
# ¿LA REGLA DE CAPACIDAD ES DE a_P, O DE CUALQUIER ESCALERA?   25 de agosto de 2026.
#
# POR QUE.  Toda la nota 2 atribuye la regla de capacidad a las POTENCIAS DE a_P.  Pero nunca hemos
# medido las potencias del OTRO elemento de Kostant, a_Q, en el mismo grupo.  Si a_Q^k cumpliera la
# misma regla, lo encontrado no seria de a_P: seria de las escaleras en general, y el encuadre de la
# nota estaria mal.  Esta gate es una FALSACION DE LO NUESTRO, y se corre antes de enviar nada.
#
# QUE ES a_Q EN EL TORO DE Sp(2m).  Se define por <alpha_i, x_Q> = 1 en toda raiz simple.  En C_m
# eso da x_i - x_{i+1} = 1 y 2 x_m = 1, luego x_Q = (2m-1, 2m-3, ..., 1)/2 = rho^v, SEMIENTERO.
# Su levantamiento al toro de Sp(2m) tiene, por tanto, orden 2h = 4m, y coordenadas
#       a_Q = ( xi^{2m-1}, xi^{2m-3}, ..., xi^1 ),      xi = zeta_{4m}
# es decir la escalera IMPAR (1,3,5,...,2m-1) evaluada en una raiz de orden t' = 4m.  (Es el
# "2 rho^v" de la lista de tres vectores con la propiedad {0,+-1}.)  Se COMPRUEBA abajo, no se supone.
#
# QUE SE MIDE
#   Q0  CONTROL: x_Q emparejado con las simples da 1 en todas; y la escalera par (a_P) NO.
#   Q1  el orden de a_Q^k, y sus autovalores.
#   Q2  el ALFABETO de a_Q^k: ¿tiene tambien la forma "d copias de mu_q menos dos puntos"?
#   Q3  LA PRUEBA DE FUEGO: ¿cumple a_Q^k la regla de capacidad de la nota?  Se prueban TODAS las
#       capacidades (cap_no_fija, cap_fija) de 0..4 y se dice cual, si alguna, acierta el 100 %.
#   Q4  y el control con a_P, que tiene que salir (d, floor(d/2)) --- si no saliera, el instrumento
#       estaria roto y Q3 no mediria nada.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python aQ_potencias.py > aQ_potencias_OUT.txt 2>&1

import cmath
import itertools
import sys
from collections import Counter
from fractions import Fraction
from math import gcd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SEP = "=" * 96
fallos = []


def ok(c, e, d=""):
    print(("   OK   " if c else "  FALLA ") + e + ("   " + d if d else ""))
    if not c:
        fallos.append(e)


def det_c(M):
    n = len(M)
    M = [r[:] for r in M]
    d = 1 + 0j
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        if abs(M[p][c]) < 1e-11:
            return 0j
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        inv = 1 / M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] * inv
            if f:
                for k in range(c, n):
                    M[r][k] -= f * M[c][k]
    return d


def dominantes(m, top):
    for lam in itertools.product(range(top + 1), repeat=m):
        if all(lam[i] >= lam[i + 1] for i in range(m - 1)):
            yield list(lam)


def h_de(expon, t, N):
    """h_n del alfabeto {zeta_t^{+-e} : e in expon}."""
    z = cmath.exp(2j * cmath.pi / t)
    h = [0j] * (N + 1)
    h[0] = 1 + 0j
    for e in expon:
        for y in (z ** (e % t), z ** ((-e) % t)):
            acc, nue = 0j, [0j] * (N + 1)
            for n in range(N + 1):
                acc = acc * y + h[n]
                nue[n] = acc
            h = nue
    return h


def chi(lam, m, h):
    H = lambda n: (h[n] if 0 <= n < len(h) else 0j)
    M = [[H(lam[i] - (i + 1) + (j + 1)) + H(lam[i] - (i + 1) - (j + 1) + 2)
          for j in range(m)] for i in range(m)]
    v = det_c(M) / 2
    if abs(v.imag) > 1e-5 or abs(v.real - round(v.real)) > 1e-5:
        return None
    return int(round(v.real))


def cls(x, q):
    return min(x % q, (-x) % q)


# ---------------------------------------------------------------- Q0
print(SEP)
print("Q0 -- CONTROL: quien es principal (Coxeter) y quien no")
print("     simples de C_m: alpha_i = e_i - e_{i+1} (i<m), alpha_m = 2 e_m")
for m in (2, 3, 4, 5):
    xQ = [Fraction(2 * (m - i) - 1, 2) for i in range(m)]      # rho^v, semientero
    xP = [Fraction(m - i) for i in range(m)]                    # la escalera par
    def pares(x):
        s = [x[i] - x[i + 1] for i in range(m - 1)] + [2 * x[m - 1]]
        return [str(v) for v in s]
    pQ, pP = pares(xQ), pares(xP)
    ok(all(v == "1" for v in pQ) and not all(v == "1" for v in pP),
       f"m={m}", f"a_Q -> {pQ}   escalera par -> {pP}")

# ---------------------------------------------------------------- Q1/Q2
print(SEP)
print("Q1/Q2 -- a_Q^k: orden y alfabeto.   t' = 4m,  exponentes k*(2j-1),  j=1..m")
print(f"{'m':>3} {'k':>3} {'tp':>4} {'d':>3} {'q':>4} {'|A|':>4} "
      f"{'A = d*mu_q - 2 ptos?':>22}   perfil del alfabeto por residuo mod q   (tp = t' = 4m)")
for m in (3, 4, 5):
    t = 4 * m
    for k in (1, 2, 3, 4):
        d = gcd(k, t)
        q = t // d
        expon = [k * (2 * j - 1) for j in range(1, m + 1)]
        z = cmath.exp(2j * cmath.pi / t)
        perfil = Counter()
        for e in expon + [-e for e in expon]:
            r = round(cmath.phase(z ** (e % t)) / (2 * cmath.pi) * q) % q
            perfil[r] += 1
        pred = Counter({r: d for r in range(q)})
        pred[0] -= 2 if d % 2 == 0 else 1
        if d % 2:
            pred[q // 2] -= 1
        pred = Counter({r: c for r, c in pred.items() if c})
        igual = perfil == pred
        print(f"{m:>3} {k:>3} {t:>4} {d:>3} {q:>4} {sum(perfil.values()):>4} "
              f"{('SI' if igual else 'NO'):>22}   {sorted(perfil.items())}")

# ---------------------------------------------------------------- Q3/Q4
print(SEP)
print("Q3/Q4 -- LA PRUEBA DE FUEGO.  ¿Que capacidades (no fija, fija) aciertan el 100 %?")
print("     Se prueban TODAS las parejas de 0..4; se imprime la lista de las que aciertan.")
print("     Para a_P la respuesta tiene que ser (d, floor(d/2)) --- eso es el control Q4.")
print()
print(f"{'elem':>5} {'m':>3} {'k':>3} {'t':>4} {'d':>3} {'q':>4} {'pesos':>6} {'nulos':>6} "
      f"{'capacidades que aciertan 100%':>34}")
for etiqueta in ("a_P", "a_Q"):
    for m in (3, 4, 5, 6):
        t = (2 * m + 2) if etiqueta == "a_P" else (4 * m)
        for k in (1, 2, 3, 4):
            d = gcd(k, t)
            q = t // d
            if q < 3:
                continue
            expon = ([k * j for j in range(1, m + 1)] if etiqueta == "a_P"
                     else [k * (2 * j - 1) for j in range(1, m + 1)])
            top = 7 if m <= 4 else 5
            h = h_de(expon, t, 2 * (top + m) + 6)
            datos = []
            for lam in dominantes(m, top):
                v = chi(lam, m, h)
                if v is None:
                    continue
                ell = [lam[i] + (m - i) for i in range(m)]
                datos.append((Counter(cls(e, q) for e in ell), v != 0))
            if not datos:
                continue
            fixed = {0} | ({q // 2} if q % 2 == 0 else set())
            aciertan = []
            for cn in range(5):
                for cf in range(5):
                    if all(all(n <= (cf if c in fixed else cn) for c, n in cnt.items()) == nz
                           for cnt, nz in datos):
                        aciertan.append((cn, cf))
            nulos = sum(1 for _, nz in datos if not nz)
            print(f"{etiqueta:>5} {m:>3} {k:>3} {t:>4} {d:>3} {q:>4} {len(datos):>6} {nulos:>6} "
                  f"{str(aciertan):>34}")
            if etiqueta == "a_P" and (d, d // 2) not in aciertan:
                fallos.append(f"Q4 control a_P m={m} k={k}")

print(SEP)
ok(not [f for f in fallos if f.startswith("Q4")],
   "CONTROL Q4: en a_P siempre acierta (d, floor(d/2))")
print()
print("LECTURA -- se escribe a la vista de la tabla, no antes.")
print(SEP)
if fallos:
    print(f"RESULTADO: {len(fallos)} FALLOS -> {sorted(set(fallos))[:6]}")
    sys.exit(1)
print("RESULTADO: el control pasa; la lectura de a_Q va en la nota.")
