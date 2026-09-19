# -*- coding: utf-8 -*-
r"""II_automata_y_conrey.py -- el automata de ceros y la identidad de la serie seno (nota II).

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

(A) AUTOMATA DE CEROS.  chi = chi_{q,+1}, q = 3 mod 4, h = h(-q), D(r) = C(r) - h.  De (D):
    S(k) + a k = b  y  k = q j + r   <=>   S(j) + (a + D(r)) j = (b - T(r) - a r)/q   (si es entero).
    Estados (a, b); inicio (0, 0) (buscamos S(m) = 0); se leen digitos desde el MENOS significativo;
    se acepta cuando j = 0 y b = 0, con ultimo digito (el mas significativo) no nulo.  BFS por
    niveles hasta DMAX digitos: 'existe un cero con <= d digitos'.  Control: q = 59 debe dar su primer
    cero 4591439 (4 digitos) y q = 11 el 11 (2 digitos).
(B) IDENTIDAD DE CONREY (Thm 4 de arXiv:2404.19647, especializada): para q = 3 mod 4,
    g_q(x) = sum chi(n) sin(2 pi n x)/n^2 = (2 pi^2 / q^{3/2}) F_q(q x),  F_q(y) = T(r) + y (h - C(r)),
    r = floor(y).  Comprobado en coma flotante (serie truncada a N terminos) en varios x.
Uso: python II_automata_y_conrey.py [DMAX]
"""
import sys, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 40


def tablas(q):
    L = [0] + [-1] * (q - 1)
    for r in range(1, q):
        L[r * r % q] = 1
    C = [0] * q; T = [0] * q; c = t = 0
    for r in range(q):
        c += L[r]; t += r * L[r]; C[r] = c; T[r] = t
    return L, C, T


def automata(q, dmax):
    L, C, T = tablas(q)
    h = -T[q - 1] // q
    D = [C[r] - h for r in range(q)]
    # estado -> menor numero leido hasta ahora (para reconstruir el primer cero)
    nivel = {(0, 0): 0}
    potencia = 1
    for d in range(1, dmax + 1):
        nuevo = {}
        mejor = None
        for (a, b), val in nivel.items():
            for r in range(q):
                num = b - T[r] - a * r
                if num % q:
                    continue
                st = (a + D[r], num // q)
                v = val + r * potencia
                if r != 0 and st[1] == 0 and (mejor is None or v < mejor):
                    mejor = v  # j = 0 aceptando: S(v) = 0 con v de d digitos
                if st not in nuevo or v < nuevo[st]:
                    nuevo[st] = v
        potencia *= q
        nivel = nuevo
        if mejor is not None:
            return d, mejor, len(nivel)
    return None, None, len(nivel)


def comprobar_directo(q, m):
    L, C, T = tablas(q)
    s = 0
    for n in range(1, m + 1):
        x, e = n, 1
        while x % q == 0:
            x //= q
        s += n * L[x % q]
    return s


for q in (11, 59, 179, 227, 347):
    d, m, estados = automata(q, DMAX if q != 347 else 16)
    extra = ""
    if m is not None and m < 10 ** 7:
        extra = "  (suma directa: S(%d) = %d)" % (m, comprobar_directo(q, m))
    print("(A) q=%3d: primer cero con %s digitos: %s ; estados en el ultimo nivel: %d%s"
          % (q, d, m, estados, extra))
    sys.stdout.flush()

# (B)
for q in (59, 179):
    L, C, T = tablas(q); h = -T[q - 1] // q
    n = np.arange(1, 400001, dtype=np.float64)
    chi = np.array([L[int(k) % q] for k in range(1, 400001)], dtype=np.float64)
    for x in (0.1, 0.23, 0.37):
        g = float(np.sum(chi * np.sin(2 * np.pi * n * x) / n ** 2))
        y = q * x; r = int(math.floor(y))
        F = T[r] + y * (h - C[r])
        print("(B) q=%d x=%.2f: g=%.8f  (2pi^2/q^1.5) F=%.8f" % (q, x, g, 2 * math.pi ** 2 / q ** 1.5 * F))
