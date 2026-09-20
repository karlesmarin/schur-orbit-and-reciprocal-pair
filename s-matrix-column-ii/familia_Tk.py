# -*- coding: utf-8 -*-
r"""familia_Tk.py -- juzgar con codigo nuestro la familia explicita de la revision.

La revision del 20-sep propone sustituir nuestra apelacion a la genericidad ("siempre que
M_3 M_5 M_7 != 0") por una construccion EXPLICITA.  Para p > 7 con (2|p) = 1 y m = (p-1)/2,
existe un unico k en [0, m] con

    8k(k+1) + 1 = 0  (mod p)        [equivale a (2k+1)^2 = 1/2]

y entonces  T_k = {1..k} U {-(k+1)..-m}  cumple, segun ellos,

    M_1 = 0,   M_3 = -1/128,   M_5 = 7/1536,   M_7 = -245/49152   (en F_p)

y de ahi  S_T = <3,5,7>,  delta = 3,  no Gorenstein, para TODO primo admisible.

Aqui no se concede nada: se construye T_k, se comprueba que es semisistema, se comparan los
tres momentos con las constantes EXACTAS que ellos dan, y el semigrupo se calcula con nuestro
s_exacto de la Parte I, que es el que usa la nota.  Tambien se mira cuantos primos admite la
condicion (2|p)=1 en el rango, porque una familia infinita que solo cubriera la mitad de los
primos sigue siendo infinita pero hay que decirlo.

    python familia_Tk.py [PMAX]

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import os
import sys
from fractions import Fraction

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(AQUI), "note_P3a", "release_tool"))
from P3a import es_primo, estabilizador, momento, s_exacto, semigrupo   # noqa: E402

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 300

CONST = {3: Fraction(-1, 128), 5: Fraction(7, 1536), 7: Fraction(-245, 49152)}


def en_Fp(fr, p):
    return fr.numerator * pow(fr.denominator, -1, p) % p


print("%-6s %-5s %-7s %-6s %-4s %-12s %-22s %s"
      % ("p", "k", "semi", "M1", "h", "huecos", "M3,M5,M7 == constantes", "veredicto"))

casos = fallos = 0
sin_k = []
for p in range(11, PMAX + 1):
    if not es_primo(p):
        continue
    if pow(2, (p - 1) // 2, p) != 1:          # (2|p) = 1  <=>  p = +-1 mod 8
        continue
    m = (p - 1) // 2
    ks = [k for k in range(0, m + 1) if (8 * k * (k + 1) + 1) % p == 0]
    if len(ks) != 1:
        sin_k.append((p, ks))
        continue
    k = ks[0]
    T = tuple(sorted([t % p for t in range(1, k + 1)]
                     + [(-t) % p for t in range(k + 1, m + 1)]))
    if len(set(T)) != m:
        print("%-6d %-5d COLISION DE ELEMENTOS" % (p, k))
        fallos += 1
        continue
    S = set(T)
    semi = 0 not in S and all((a in S) != ((p - a) in S) for a in range(1, (p + 1) // 2))
    M1 = momento(T, 1, p)
    h = len(estabilizador(T, p))
    try:
        sg = semigrupo(s_exacto(T, p))
    except ValueError:
        print("%-6d %-5d s_exacto no concluye" % (p, k))
        fallos += 1
        continue
    huecos = tuple(sg["huecos"])
    consts_ok = all(momento(T, j, p) == en_Fp(CONST[j], p) for j in (3, 5, 7))
    ok = semi and M1 == 0 and h == 1 and huecos == (1, 2, 4) and consts_ok
    casos += 1
    fallos += 0 if ok else 1
    print("%-6d %-5d %-7s %-6d %-4d %-12s %-22s %s"
          % (p, k, "si" if semi else "NO", M1, h, str(huecos), "si" if consts_ok else "NO",
             "ok" if ok else "FALLA"))
    sys.stdout.flush()

print("")
print("primos admisibles probados: %d ; fallos: %d" % (casos, fallos))
if sin_k:
    print("primos sin k unico en [0,m]: %s" % sin_k[:6])
print("NOTA DE ALCANCE: la construccion pide (2|p) = 1, o sea p = +-1 mod 8, que es la mitad de")
print("los primos por densidad de Dirichlet; la familia es infinita, pero no cubre todo primo.")
print("VEREDICTO:", "la construccion de la revision se sostiene; convierte la familia generica "
      "en una proposicion demostrable" if not fallos and casos else "ALGO FALLA")
