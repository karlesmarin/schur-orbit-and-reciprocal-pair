# -*- coding: utf-8 -*-
r"""frontera_efectiva.py -- la curva de Bring generalizada para n = 6: ¿tiene tambien nombre?

DONDE NADIE MIRA.  Para todo n, anular M_1 = ... = M_{n-2} equivale a anular e_1,...,e_{n-2}, y el
polinomio queda X^n - e_{n-1} X + e_n: un TRINOMIO.  Asi que las columnas de defecto maximo son,
para cada n, los trinomios X^n + aX + b totalmente escindidos, y el lugar correspondiente es una
curva C_n en P^{n-1}.  Para n = 5 es la curva de Bring y su jacobiana es E^4 con E = 50.a3
(verificado en 51 primos, extremo_bring.py).  Nadie ha mirado n = 6.

QUE HACE ESTE GUION.  Cuenta, para cada primo, los pares (a,b) con a,b != 0 tales que X^6 + aX + b
se escinde con raices distintas, reconstruye #C_6(F_p) sumando los tres tipos de frontera con el
mismo argumento que en n = 5, y saca la "traza" t_p = p + 1 - #C_6(F_p).  Luego busca estructura:
tamano frente a la cota de Hasse-Weil con el genero que da la formula de adjuncion, divisibilidad,
y si t_p es multiplo de trazas de curvas pequenas.

El genero por adjuncion (su formula, seccion 6): g_n = 1 + [(n-2)(n-3) - 4] (n-2)! / 4.
Para n = 6: g = 1 + (12 - 4)*24/4 = 49.

    python frontera_efectiva.py [PMAX]

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import math
import os
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(AQUI), "note_P3a", "release_tool"))
from P3a import es_primo                                              # noqa: E402
sys.path.insert(0, AQUI)
from extremo_bring import escinde_distinto                     # noqa: E402

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 120
N = 6
G = 1 + ((N - 2) * (N - 3) - 4) * math.factorial(N - 2) // 4


def R_escinde(p, n):
    """¿se escinde R_n(X) = (X^n - nX + n-1)/(X-1)^2 sobre F_p?  (raices distintas)"""
    # coeficientes de X^n - nX + (n-1) divididos por (X-1)^2, por division sintetica
    c = [0] * (n + 1)
    c[n], c[1], c[0] = 1, -n, n - 1
    for _ in range(2):
        q, r = [0] * (len(c) - 1), 0
        for i in range(len(c) - 1, -1, -1):
            r = r + c[i]
            if i:
                q[i - 1] = r
            r = r * 1
        # division por (X - 1): coeficientes de q
        c = q
    grado = len(c) - 1
    raices = [x for x in range(p) if sum(ci * pow(x, i, p) for i, ci in enumerate(c)) % p == 0]
    return len(set(raices)) == grado


print("n = %d, genero por adjuncion g = %d, cota de Hasse-Weil 2g sqrt(p) = %.1f sqrt(p)" % (N, G, 2 * G))
print("%-6s %-10s %-12s %-10s %-12s %s" % ("p", "N_orb", "#C_6(F_p)", "t_p", "t_p/sqrt(p)", "|t|<=2g sqrt(p)"))
datos = []
for p in range(7, PMAX + 1):
    if not es_primo(p):
        continue
    pares = sum(1 for a in range(1, p) for b in range(1, p) if escinde_distinto(a, b, p, N))
    N_orb = pares // (p - 1)
    epsR = 1 if R_escinde(p, N) else 0
    eps5 = 1 if (p - 1) % (N - 1) == 0 else 0
    eps6 = 1 if (p - 1) % N == 0 else 0
    # #C_n = n! * [ N_orb + epsR/2 + eps_{n-1}/(n-1) + eps_n/n ]
    C = math.factorial(N) * N_orb + math.factorial(N) // 2 * epsR + \
        math.factorial(N) // (N - 1) * eps5 + math.factorial(N) // N * eps6
    t = p + 1 - C
    datos.append((p, t))
    print("%-6d %-10d %-12d %-10d %-12.3f %s"
          % (p, N_orb, C, t, t / math.sqrt(p), "si" if abs(t) <= 2 * G * math.sqrt(p) else "NO"))
    sys.stdout.flush()

print("")
print("=" * 96)
print("LA FRONTERA EFECTIVA: hasta donde sirve la geometria")
print("=" * 96)
print("El termino principal del recuento de columnas de defecto maximo es p/n!, y la cota de")
print("Hasse-Weil es 2 g_n sqrt(p).  El recuento solo es informativo cuando p > (2 g_n)^2.")
print("")
print("%-4s %-12s %-10s %-16s %s" % ("n", "genero g_n", "2 g_n", "densidad 1/n!", "p necesario"))
for n in range(5, 10):
    g = 1 + ((n - 2) * (n - 3) - 4) * math.factorial(n - 2) // 4
    print("%-4d %-12d %-10d 1/%-14d %d" % (n, g, 2 * g, math.factorial(n), (2 * g) ** 2))
print("")
print("Es decir: n = 5 es el ULTIMO caso en que la curva sirve (p > 64, y de hecho el primer")
print("primo con columna de defecto maximo es 233).  Para n = 6 haria falta p > 9604; para n = 7,")
print("p > 925444.  El genero crece como n!/4, asi que la cota se vuelve vacia mas rapido de lo que")
print("crece el termino principal.  Y la cota de Weil para sumas de caracteres, H_r(n,p) =")
print("binom((r-1)sqrt(p)+n-1, n), degenera igual de rapido.  La geometria NO da la asintotica")
print("para n general: es una herramienta de n = 5.")
print("")
print("estructura de t_p:")
print("   todos pares       :", all(t % 2 == 0 for _, t in datos))
print("   divisibles por 4  :", all(t % 4 == 0 for _, t in datos))
print("   divisibles por 6  :", all(t % 6 == 0 for _, t in datos))
print("   maximo |t|/sqrt(p):", max(abs(t) / math.sqrt(p) for _, t in datos))
print("   secuencia t_p     :", [t for _, t in datos])
