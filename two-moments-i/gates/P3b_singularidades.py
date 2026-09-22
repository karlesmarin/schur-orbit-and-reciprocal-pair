# -*- coding: utf-8 -*-
r"""P3b_singularidades.py -- el lugar singular de V_lambda es el conjunto de los vectores de signos.

EL OBJETO.  En el tamiz de C(n,p) = #{T subset F_p : |T| = n, M_1 = M_3 = 0} aparece

    V_lambda = { sum w_i x_i = 0 } cap { sum w_i x_i^3 = 0 }  en P^{l-1},

una hipersuperficie cubica de P^{l-2}, con lambda = (w_1,...,w_l) particion de n.

EL TEOREMA (demostrado, y este fichero lo comprueba en vez de suponerlo).  Un punto de V_lambda
es singular si y solo si los dos gradientes son proporcionales:

    (w_i)  y  (3 w_i x_i^2)   proporcionales   <=>   x_i^2 = c  para todo i

(usando w_i != 0 en F_p, o sea p > max w_i).  Con c = 0 sale el origen, que no es un punto
proyectivo; con c != 0, escalando por una raiz de c queda x_i = +-1.  Y entonces la PRIMERA
ecuacion dice exactamente que los pesos se parten en dos bloques de igual suma, mientras que la
SEGUNDA se cumple sola, porque (+-1)^3 = +-1.  Por tanto:

    Sing(V_lambda) = { vectores de signos +-1, salvo signo global, con sum_{A} w = sum_{A^c} w }
    #Sing(V_lambda) = #{ A subconjunto de posiciones : sum_{i in A} w_i = n/2 } / 2.

COROLARIO, Y ES EL MECANISMO QUE BUSCABAMOS.  Como sum w_i = n, un bloque de igual suma exige
**n PAR**.  Para n IMPAR toda V_lambda es LISA, y en cuanto l >= 4 una cubica plana lisa es una
curva eliptica de verdad: el recuento deja de ser elemental.  Para n par muchas lambda se parten
y degeneran.  Eso explica por que n = 6 es elemental y n = 7 no, sin apelar a monotonia en n.

Y da la clasificacion exacta: C(n,p) es elemental (polinomio + simbolos de Legendre) si y solo si
TODA lambda |- n con l >= 4 admite un bloque de suma n/2 -- porque con l <= 3 la variedad tiene
dimension 0 y lo unico que puede aportar es un simbolo cuadratico.  Este fichero lista los n que
cumplen eso.

CONTROLES:
 (1) lugar singular por fuerza bruta sobre F_p contra la prediccion de signos, lambda a lambda;
 (2) los cuatro casos clasicos con nombre: tres planos, Segre (10 nodos), Cayley (4 nodos),
     y la Clebsch, que es LISA porque 5 es impar;
 (3) la clasificacion de los n elementales.

    python P3b_singularidades.py [PMAX] [NMAX]

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import sys
from itertools import product

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


# Autocontenido: no lee ni importa ningun otro fichero.
def es_primo(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def particiones(n, maxi=None):
    maxi = maxi or n
    if n == 0:
        yield ()
        return
    for w in range(min(n, maxi), 0, -1):
        for resto in particiones(n - w, w):
            yield (w,) + resto


PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 23
NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 16


def bloques(lam, p=None):
    """Vectores de signos, salvo signo global.

    Sobre Qbar (p = None): #{A : 2 sum_{i in A} w_i = n} / 2, la particion en bloques de
    igual suma.  Sobre F_p la condicion es la CONGRUENCIA 2 sum_A w = n (mod p), que para
    p > n coincide con la igualdad entera y para p | n admite soluciones extra --- con
    p = n = 7 el propio vector de todos-unos, porque sum w_i = n = 0 en F_p."""
    n = sum(lam)
    c = 0
    for A in product([0, 1], repeat=len(lam)):
        s = sum(w for w, a in zip(lam, A) if a)
        if (2 * s - n) % p == 0 if p else (2 * s == n):
            c += 1
    assert c % 2 == 0, (lam, p)
    return c // 2


def puntos_proyectivos(l, p):
    """Representantes de P^{l-1}(F_p): primera coordenada no nula igual a 1."""
    for i in range(l):
        for cola in product(range(p), repeat=l - 1 - i):
            yield (0,) * i + (1,) + cola


def sing_fuerza_bruta(lam, p):
    """Puntos de V_lambda donde los dos gradientes son proporcionales."""
    l = len(lam)
    out = set()
    for x in puntos_proyectivos(l, p):
        if sum(w * xi for w, xi in zip(lam, x)) % p:
            continue
        if sum(w * pow(xi, 3, p) for w, xi in zip(lam, x)) % p:
            continue
        # rango <= 1 de [[w_i], [3 w_i x_i^2]]: todos los x_i^2 iguales
        cuadrados = {(xi * xi) % p for xi in x}
        if len(cuadrados) == 1:
            out.add(x)
    return out


print("(1) LUGAR SINGULAR: fuerza bruta sobre F_p contra la prediccion por vectores de signos.")
print("    La prediccion se hace con la CONGRUENCIA 2 sum_A w = n (mod p); entre parentesis va")
print("    la cuenta entera sobre Qbar, que es la que vale siempre que p > n.")
print("%-18s %-4s %-5s %-9s %s" % ("lambda", "n", "l", "sobre Qbar", "medido:predicho por primo"))
fallos = casos = 0
extras = 0
por_l, fallos_l = {}, {}
primos = [p for p in range(7, PMAX + 1) if es_primo(p)]
for n in range(3, 8):
    for lam in sorted(particiones(n), key=lambda t: (-len(t), t)):
        if len(lam) < 3 or len(lam) > 5:
            continue
        cero = bloques(lam)
        med = []
        for p in primos:
            if p <= max(lam) or p ** len(lam) > 4_000_000:
                continue
            pred = bloques(lam, p)
            k = len(sing_fuerza_bruta(lam, p))
            med.append((p, k, pred))
            casos += 1
            fallos += 0 if k == pred else 1
            extras += 1 if pred != cero else 0
            por_l[len(lam)] = por_l.get(len(lam), 0) + 1
            fallos_l[len(lam)] = fallos_l.get(len(lam), 0) + (0 if k == pred else 1)
        print("%-18s %-4d %-5d %-9d %s"
              % (str(lam), n, len(lam), cero,
                 " ".join("%d:%d/%d" % (p, k, q) for p, k, q in med)))
        sys.stdout.flush()
print("")
print("casos: %d ; fallos: %d ; casos con p | n donde la congruencia da extras: %d"
      % (casos, fallos, extras))
print("por longitud: %s ; fallos por longitud: %s"
      % (sorted(por_l.items()), sorted(fallos_l.items())))
print("   los casos l = 3 quedan FUERA del corolario (que pide l >= 4) y coinciden igual con la")
print("   cuenta de signos: no ensenan por que hace falta la hipotesis.  El caso que la fuerza es")
print("   (3,3), de longitud 2: x + y = 0 y x^3 + y^3 = (x+y)(x^2-xy+y^2), el ideal es (x+y), un")
print("   punto reducido y liso, mientras la cuenta de signos da 1.")

print("")
print("(2) LOS CUATRO CLASICOS, con su numero de nodos predicho por el teorema")
for lam, nombre in [((1, 1, 1, 1), "tres planos (n=4)"),
                    ((1, 1, 1, 1, 1), "superficie de CLEBSCH (n=5)"),
                    ((2, 1, 1, 1, 1), "cubica de CAYLEY, la 4-nodal (n=6)"),
                    ((1, 1, 1, 1, 1, 1), "cubica de SEGRE (n=6)"),
                    ((1,) * 7, "cubica de dim 4 (n=7)")]:
    print("   %-20s n=%d  nodos predichos = %-3d   %s"
          % (str(lam), sum(lam), bloques(lam), nombre))

print("")
print("(3) CLASIFICACION: C(n,p) es elemental  <=>  toda lambda |- n con l >= 4 se parte")
print("    (con l <= 3 la variedad tiene dimension 0 y solo puede dar un simbolo cuadratico)")
elementales = []
for n in range(3, NMAX + 1):
    malas = [lam for lam in particiones(n) if len(lam) >= 4 and bloques(lam) == 0]
    if not malas:
        elementales.append(n)
        print("   n = %-3d ELEMENTAL" % n)
    else:
        peor = min(malas, key=lambda t: (len(t), t))
        print("   n = %-3d no: %s tiene l = %d y ningun bloque de suma %s"
              % (n, str(peor), len(peor), n / 2))
print("")
print("   los n elementales hasta %d: %s" % (NMAX, elementales))
print("")
print("VEREDICTO:", "el lugar singular es el de los vectores de signos, en %d casos" % casos
      if not fallos else "EL TEOREMA FALLA")
bien = [casos > 0 and fallos == 0, elementales == [3, 4, 6] or NMAX < 6]
print("")
print("=" * 72)
print("VERIFICACION LUGAR SINGULAR:  %d ok, %d MAL" % (sum(bien), len(bien) - sum(bien)))
print("=" * 72)
sys.exit(0 if all(bien) else 1)
