# -*- coding: utf-8 -*-
r"""cascada_newton.py -- por que en n = 4 no hay delta = 2, y la recurrencia que lo explica.

EL HECHO MEDIDO (censo_delta.py).  En n = 4 el reparto de delta es {0, 1, 3}: el valor 2 NO
aparece para ningun p <= 31.  En n = 5, 6, 7 si aparece.

LA EXPLICACION.  Los momentos de un conjunto T con |T| = n cumplen la recurrencia de Newton
    M_r = e_1 M_{r-1} - e_2 M_{r-2} + ... + (-1)^{n-1} e_n M_{r-n}   para r > n,
porque e_j = 0 para j > n.  Con T centrado (e_1 = 0) y ademas e_2 = 0 (es decir M_2 = 0), para
n = 4 se obtiene
    M_5 = e_1 M_4 - e_2 M_3 + e_3 M_2 - e_4 M_1 = 0,
asi que M_2 = 0 arrastra M_5 = 0: el semigrupo pierde a la vez el 2 y el 5 y queda <3,4>, con
huecos {1,2,5} y delta = 3.  Nunca delta = 2.  El mismo argumento dice que la anulacion de momentos
no es independiente: los momentos son una sucesion recurrente lineal sobre F_p, y S_T es el
soporte de esa sucesion.

ESTE GUION lo comprueba en vez de creerlo: verifica la recurrencia en todas las orbitas, confirma
la implicacion M_2 = 0 => M_5 = 0 en n = 4, y busca delta = 2 en n = 4 (no debe existir).

    python cascada_newton.py [PMAX]

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import os
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(AQUI), "note_P3a", "release_tool"))
from P3a import centrados, es_primo, momento, semigrupo, s_exacto        # noqa: E402

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 37


def elementales(T, p):
    """e_0 .. e_n de T en F_p, por el polinomio prod (X - t)."""
    coef = [1]
    for t in T:
        coef = [(c - (coef[i - 1] if i else 0) * t) % p for i, c in enumerate(coef + [0])]
    # coef[i] = coeficiente de X^{n-i} = (-1)^i e_i
    return [((-1) ** i * coef[i]) % p for i in range(len(coef))]


mal_rec, mal_casc, delta2_n4, casos = 0, 0, [], 0
for p in range(5, PMAX + 1):
    if not es_primo(p):
        continue
    for n in range(3, min(7, p - 2) + 1):
        for T in centrados(p, n):
            casos += 1
            e = elementales(T, p)
            M = [n % p] + [momento(T, r, p) % p for r in range(1, 3 * n + 2)]
            for r in range(n + 1, 3 * n + 2):                 # la recurrencia, r > n
                pred = sum((-1) ** (i - 1) * e[i] * M[r - i] for i in range(1, n + 1)) % p
                if pred != M[r]:
                    mal_rec += 1
            if n == 4:
                if M[2] == 0 and M[5] != 0:                   # la cascada
                    mal_casc += 1
                if semigrupo(s_exacto(T, p))["genero"] == 2:
                    delta2_n4.append((p, T))

print("orbitas examinadas: %d  (p <= %d, 3 <= n <= 7)" % (casos, PMAX))
print("fallos de la recurrencia de Newton      : %d" % mal_rec)
print("fallos de la cascada M_2=0 => M_5=0 (n=4): %d" % mal_casc)
print("orbitas con n = 4 y delta = 2            : %d %s" % (len(delta2_n4), delta2_n4[:5]))
print("")
print("VEREDICTO:", "la explicacion SOBREVIVE" if not (mal_rec or mal_casc or delta2_n4)
      else "ALGO FALLA, leer arriba")
