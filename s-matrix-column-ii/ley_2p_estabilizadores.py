# -*- coding: utf-8 -*-
r"""ley_2p_estabilizadores.py -- la correccion que le faltaba a la tabla de F_8.

EL FALLO QUE REPARA (senalado en la revision del 20-sep y confirmado aqui).  ley_2p_largo.py
calcula, y su propia cabecera lo dice,

    F = (C12 + C13 - C123) / (C(p,n)/p),

que es la proporcion de espectros centrados en el LOCUS DE MOMENTOS  M_2 M_3 = 0.  Pero el
teorema define F_n(p) como la proporcion con delta >= 2, y el Corolario 3 solo identifica las
dos cosas CUANDO h = 1.  El cuarto paso de la demostracion existe justamente para controlar esa
diferencia, y la tabla la omitia.  O sea la tabla y la figura 5 median otra variable.

LA CORRECCION ES EXACTA Y ENUMERABLE, no hace falta estimarla.  Clasificacion completa de los
espectros centrados de tamano n = 8 con estabilizador no trivial:

    H = Stab(T) de orden h > 1.  T menos el cero es union de H-clases, todas de tamano h.
      * si 0 no esta en T:  h | 8, o sea h en {2,4,8}.  En los tres casos -1 esta en H, luego
        T es una UNION DE CUATRO PAREJAS {t,-t}.  Son C((p-1)/2, 4) conjuntos.
      * si 0 esta en T:  h | 7, o sea h = 7, y hace falta 7 | p-1.  Son (p-1)/7 conjuntos.

No hay mas casos.  Para p = 31 eso da C(15,4) = 1365 y ningun caso de orden 7, que es
exactamente el numero que la revision reporta como sobrante.

Para cada uno se decide delta de verdad: 1 pertenece a S_T si y solo si M_h != 0, y si M_h != 0
entonces S_T = N_0 y delta = 0 --- eso resuelve el 99% con UN momento.  Solo cuando M_h = 0 se
calcula S_T entero con s_exacto de la Parte I.

    python ley_2p_estabilizadores.py [PMAX]

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import os
import sys
from itertools import combinations
from math import comb

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(AQUI), "note_P3a", "release_tool"))
from P3a import es_primo, estabilizador, momento, s_exacto, semigrupo  # noqa: E402

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 97
N = 8


def raiz_primitiva(p):
    fac = [q for q in range(2, p) if (p - 1) % q == 0 and es_primo(q)]
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in fac):
            return g
    raise AssertionError(p)


def divisores(m):
    return [d for d in range(1, m + 1) if m % d == 0]


def corrige(p, verbose=False):
    """(numero de T con h>1, cuantos de ellos tienen delta>=2, control de la clasificacion)."""
    m = (p - 1) // 2
    g = raiz_primitiva(p)
    pot = [pow(g, i, p) for i in range(p - 1)]
    divs = divisores(p - 1)

    total = con_delta2 = 0
    reparto = {}

    # ---- caso A: uniones de cuatro parejas {t,-t}
    for S in combinations(range(m), 4):
        T = tuple(sorted([pot[i] for i in S] + [pot[i + m] for i in S]))
        A = set(S) | {i + m for i in S}
        h = max(d for d in divs if all(((a + (p - 1) // d) % (p - 1)) in A for a in A))
        total += 1
        if momento(T, h, p):                 # 1 pertenece a S_T  =>  S_T = N_0
            d = 0
        else:
            d = semigrupo(s_exacto(T, p))["genero"]
        reparto[d] = reparto.get(d, 0) + 1
        if d >= 2:
            con_delta2 += 1

    # ---- caso B: clase del subgrupo de orden 7, mas el cero
    if (p - 1) % 7 == 0:
        H = [pow(g, (p - 1) // 7 * i, p) for i in range(7)]
        # UN representante por clase: g^j con j = 0 .. (p-1)/7 - 1.  Recorrer todo c de 1 a p-1
        # da cada clase SIETE veces, y el control lo canto (1029 frente a 1005 en p = 29).
        for j in range((p - 1) // 7):
            c = pow(g, j, p)
            T = tuple(sorted({0} | {c * x % p for x in H}))
            assert len(T) == 8 and sum(T) % p == 0, (p, T)
            total += 1
            d = semigrupo(s_exacto(T, p))["genero"]
            reparto[d] = reparto.get(d, 0) + 1
            if d >= 2:
                con_delta2 += 1

    # ---- control: la clasificacion dice cuantos deberia haber
    esperado = comb(m, 4) + ((p - 1) // 7 if (p - 1) % 7 == 0 else 0)
    return total, con_delta2, esperado, reparto


print("n = %d.  U = C12 + C13 - C123 cuenta el LOCUS; F_8 cuenta delta >= 2." % N)
print("Todo T centrado con h > 1 esta en el locus (h par => todos los momentos impares se")
print("anulan => M_3 = 0), asi que la correccion es siempre NEGATIVA o nula.")
print("")
print("%-5s %-10s %-9s %-9s %-11s %-13s %s"
      % ("p", "h>1", "de esos d>=2", "control", "correccion", "(F-2/p)p^2+1", "reparto de delta"))
for p in range(29, PMAX + 1):
    if not es_primo(p):
        continue
    tot, d2, esp, rep = corrige(p)
    assert tot == esp, (p, tot, esp)
    print("%-5d %-10d %-9d %-9s %-11d %-13s %s"
          % (p, tot, d2, "ok" if tot == esp else "MAL", d2 - tot, "(ver abajo)",
             dict(sorted(rep.items()))))
    sys.stdout.flush()
print("")
print("La correccion es (numero con delta>=2) - (numero en el locus) = d2 - h>1, siempre <= 0.")
print("Numerador correcto de F_8(p) = U + correccion.")

# ---------------------------------------------------------------------------------------------
# LA TABLA CORREGIDA, juntando esto con el censo de tres momentos ya archivado.  No se repite el
# censo --- cuesta O(n p^4) --- sino que se lee su salida y se le aplica la correccion exacta.
print("")
print("=" * 96)
print("LA TABLA DE LA NOTA, ANTES Y DESPUES.  La publicada media el locus; la buena mide F_8.")
print("=" * 96)
import io                                                               # noqa: E402
import os as _os                                                        # noqa: E402

SAL = _os.path.join(AQUI, "ley_2p_largo_OUT.txt")
U = {}
for linea in io.open(SAL, encoding="utf-8", errors="replace"):
    f = linea.split()
    if len(f) == 7 and f[0].isdigit():
        U[int(f[0])] = (int(f[1]), int(f[2]), int(f[3]))

print("%-5s %-10s %-11s %-11s %-14s %-14s %s"
      % ("p", "U (locus)", "correccion", "numerador", "publicado", "corregido", "signo"))
for p in sorted(U):
    if p > PMAX:
        continue
    c12, c13, c123 = U[p]
    u = c12 + c13 - c123
    tot, d2, esp, rep = corrige(p)
    num = u + (d2 - tot)
    cent = comb(p, N) // p
    viejo = (u / cent - 2 / p) * p * p + 1
    nuevo = (num / cent - 2 / p) * p * p + 1
    print("%-5d %-10d %-11d %-11d %-14.4f %-14.4f %s"
          % (p, u, d2 - tot, num, viejo, nuevo, "cambia de signo" if viejo * nuevo < 0 else ""))
    sys.stdout.flush()
print("")
print("El residuo corregido se acerca a cero POR ABAJO, no por arriba: la nota lo contaba al reves.")
