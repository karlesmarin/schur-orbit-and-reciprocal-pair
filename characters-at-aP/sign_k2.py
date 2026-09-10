# -*- coding: utf-8 -*-
# EL SIGNO EN k=2.   8 de septiembre de 2026.
#
# POR QUE.  El articulo da el VALOR ABSOLUTO en k=2 --- |chi| = prod_c D_c --- y dice, literalmente,
# que el signo no lo ha determinado.  Es el primero de los tres problemas que deja abiertos.  Esta
# gate intenta cerrarlo, y esta escrita para que FALLE si el resultado es un ajuste y no una ley.
#
# EL MONTAJE.  En Sp(2m):  t = 2m+2,  k = 2,  d = 2,  q = m+1.  El alfabeto de a_P^2 es
# 2*mu_q menos dos copias de 1, luego  H(z) = (1-z)^2 / (1-z^q)^2  y los h_n son ENTEROS: el
# caracter se calcula exacto, sin flotantes, y su signo es un dato duro.
#
# EL METODO, y es lo unico interesante de la gate.  Se calculan quince estadisticos enteros de
# lambda --- cuantas clases usan el signo mas, cuantas el menos, la paridad de los cocientes, las
# inversiones de la palabra de clases, etc. --- y se busca un subconjunto S tal que
#        signo(chi) = (-1)^{suma de los estadisticos de S}
# resolviendo un sistema LINEAL SOBRE GF(2).  Tres salvaguardas contra el autoengano:
#   T1  SENUELO: una columna de paridades ALEATORIAS entra en el sistema.  Si la solucion la usa,
#       lo que hay es sobreajuste y la gate lo dice.
#   T2  ENTRENAMIENTO Y PRUEBA SEPARADOS: la ley se ajusta con m <= 6 y se VALIDA en m = 7, 8, 9,
#       que el ajuste no ha visto.  Una ley que no cruce esa frontera no es una ley.
#   T3  CONTROL DEL INSTRUMENTO: |chi| calculado tiene que coincidir con prod_c D_c del articulo en
#       todos los supervivientes.  Si no coincide, no se mide el signo de nada.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python sign_k2.py > sign_k2_OUT.txt 2>&1

import random
import sys
from collections import Counter
from itertools import combinations

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SEP = "=" * 96
random.seed(20260908)


def dominants(m, top):
    def rec(k, hi):
        if k == 0:
            yield ()
            return
        for v in range(hi, -1, -1):
            for rest in rec(k - 1, v):
                yield (v,) + rest
    return rec(m, top)


def h_series(q, N):
    """H(z) = (1-z)^2 / (1-z^q)^2, coeficientes enteros."""
    base = [0] * (N + 1)
    k = 0
    while q * k <= N:
        base[q * k] = k + 1
        k += 1
    h = [0] * (N + 1)
    for n in range(N + 1):
        v = base[n]
        if n >= 1:
            v -= 2 * base[n - 1]
        if n >= 2:
            v += base[n - 2]
        h[n] = v
    return h


def det_int(M):
    """Bareiss: determinante exacto de una matriz entera."""
    n = len(M)
    A = [row[:] for row in M]
    sign = 1
    prev = 1
    for c in range(n):
        if A[c][c] == 0:
            piv = next((r for r in range(c + 1, n) if A[r][c] != 0), None)
            if piv is None:
                return 0
            A[c], A[piv] = A[piv], A[c]
            sign = -sign
        for r in range(c + 1, n):
            for kk in range(c + 1, n):
                A[r][kk] = (A[r][kk] * A[c][c] - A[r][c] * A[c][kk]) // prev
            A[r][c] = 0
        prev = A[c][c]
    return sign * A[n - 1][n - 1]


def chi(lam, m, h):
    H = lambda n: (h[n] if 0 <= n < len(h) else 0)
    M = [[H(lam[i] - (i + 1) + (j + 1)) + H(lam[i] - (i + 1) - (j + 1) + 2)
          for j in range(m)] for i in range(m)]
    return det_int(M) // 2


def cls(x, q):
    return min(x % q, (-x) % q)


def analiza(lam, m):
    """devuelve (superviviente?, |chi| predicho, estadisticos) para k=2."""
    q = m + 1
    ell = [lam[i] + (m - i) for i in range(m)]
    fijas = {0} | ({q // 2} if q % 2 == 0 else set())
    porclase = {}
    for i, e in enumerate(ell):
        porclase.setdefault(cls(e, q), []).append(e)
    for c, filas in porclase.items():
        if len(filas) > (1 if c in fijas else 2):
            return False, None, None, None
    pred = 1
    n_mas = n_menos = n_fija1 = 0
    for c, filas in porclase.items():
        if c in fijas:
            if len(filas) == 1:
                pred *= (2 * filas[0]) // q
                n_fija1 += 1
        elif len(filas) == 2:
            a, b = filas
            if (a - b) % q == 0:
                pred *= abs(a - b) // q
                n_menos += 1
            else:
                pred *= (a + b) // q
                n_mas += 1
    clases = [cls(e, q) for e in ell]
    inv = sum(1 for i in range(m) for j in range(i + 1, m) if clases[i] < clases[j])
    est = {
        "s1_mas": n_mas,
        "s2_menos": n_menos,
        "s3_fija1": n_fija1,
        "s4_cocientes": sum(e // q for e in ell),
        "s5_dobladas": sum(1 for e in ell if e % q > q / 2),
        "s6_|lambda|": sum(lam),
        "s7_i*lambda": sum((i + 1) * lam[i] for i in range(m)),
        "s8_inv": inv,
        "s9_m": m,
        "s10_q": q,
        "s11_libre": min(set(range(q // 2 + 1)) - set(clases)) if set(range(q // 2 + 1)) - set(clases) else 0,
        "s12_(lam,rho)": sum(lam[i] * (m - i) for i in range(m)),
        "s13_sum_ell": sum(ell),
        "s14_dobles": sum(1 for c, f in porclase.items() if len(f) == 2),
        "s15_SENUELO": random.randint(0, 1),
    }
    ocup = tuple(sorted((c, len(f)) for c, f in porclase.items()))
    return True, pred, est, ocup


NOMBRES = None
print(SEP)
print("T3 -- CONTROL DEL INSTRUMENTO y recogida de datos")
DATOS = {}
for m in range(2, 10):
    q = m + 1
    top = 10 if m <= 5 else (7 if m <= 7 else 5)
    h = h_series(q, 3 * (top + m) + 6)
    filas = []
    malos = 0
    for lam in dominants(m, top):
        vivo, pred, est, patron = analiza(lam, m)
        v = chi(lam, m, h)
        if not vivo:
            if v != 0:
                malos += 1
            continue
        if abs(v) != pred:
            malos += 1
            continue
        if NOMBRES is None:
            NOMBRES = sorted(est)
        filas.append(([est[k] % 2 for k in NOMBRES], 1 if v < 0 else 0, patron, lam))
    DATOS[m] = filas
    print(f"   m={m:>2}  q={q:>2}  supervivientes {len(filas):>7}   discrepancias |chi| vs prod D_c: {malos}")
    if malos:
        print("  FALLA el control del instrumento: no se mide el signo de nada.")
        sys.exit(1)

print(SEP)
print("T4 -- ¿DEPENDE EL SIGNO SOLO DE LA FORMA?   patron = que clases llevan 0, 1 o 2 filas")
grupos = {}
for m in DATOS:
    for f, t, patron, lam in DATOS[m]:
        grupos.setdefault((m, patron), []).append((t, lam))
consts = sum(1 for g in grupos.values() if len({t for t, _ in g}) == 1)
print(f"   {len(grupos)} patrones distintos; el signo es constante en {consts} de ellos")
if consts < len(grupos):
    ej = next((k, v) for k, v in grupos.items() if len({t for t, _ in v}) > 1)
    (mm, pat), casos = ej
    a = next(l for t, l in casos if t == 0)
    b = next(l for t, l in casos if t == 1)
    print(f"   CONTRAEJEMPLO minimo:  m={mm}, patron {pat}")
    print(f"      lambda={a} da signo +   y   lambda={b} da signo -")
    print("   => el signo NO es funcion del patron de ocupacion: mira los numeros, no la forma.")
else:
    print("   => el signo ES funcion del patron de ocupacion.  Eso lo cierra en forma cerrada.")

print(SEP)
print("T1/T2 -- LEY DEL SIGNO: se ajusta con m<=6 y se valida en m=7,8,9")


def resolver(filas, ncols):
    """sistema lineal sobre GF(2): busca x con  fila . x + const = objetivo."""
    A = [r[0][:] + [1, r[1]] for r in filas]       # columnas + constante + objetivo
    n = ncols + 1
    piv = []
    r = 0
    for c in range(n):
        p = next((i for i in range(r, len(A)) if A[i][c]), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        for i in range(len(A)):
            if i != r and A[i][c]:
                A[i] = [(a ^ b) for a, b in zip(A[i], A[r])]
        piv.append(c)
        r += 1
        if r == len(A):
            break
    for i in range(r, len(A)):
        if A[i][n] and not any(A[i][:n]):
            return None                              # inconsistente
    x = [0] * n
    for i, c in enumerate(piv):
        x[c] = A[i][n]
    return x


TRAIN = [f for m in (2, 3, 4, 5, 6) for f in DATOS[m]]
TEST = [(m, DATOS[m]) for m in (7, 8, 9)]
sol = resolver(TRAIN, len(NOMBRES))
if sol is None:
    print("   NO hay ninguna combinacion lineal de estas paridades que de el signo.")
    print("   RESULTADO: el signo NO es una paridad de estos quince estadisticos.  Queda abierto,")
    print("              y ahora se sabe que no esta aqui.")
    sys.exit(0)

usados = [NOMBRES[i] for i in range(len(NOMBRES)) if sol[i]]
print(f"   ley candidata:  signo = (-1)^( {' + '.join(usados) if usados else '0'}"
      f"{' + 1' if sol[-1] else ''} )")
if "s15_SENUELO" in usados:
    print("  FALLA T1: la ley usa la columna SENUELO -> es sobreajuste, no ley.")
    sys.exit(1)
print("   OK  T1: la ley NO usa el senuelo aleatorio.")


def comprueba(filas, x):
    mal = 0
    for f, t, _p, _l in filas:
        s = x[-1]
        for i, v in enumerate(f):
            if x[i] and v:
                s ^= 1
            elif x[i]:
                pass
        s = (sum(x[i] * f[i] for i in range(len(f))) + x[-1]) % 2
        if s != t:
            mal += 1
    return mal


mal_train = comprueba(TRAIN, sol)
print(f"   entrenamiento (m=2..6): {len(TRAIN)} casos, {mal_train} fallos")
total_mal = 0
for m, filas in TEST:
    mal = comprueba(filas, sol)
    total_mal += mal
    print(f"   VALIDACION m={m}: {len(filas):>7} casos, {mal:>5} fallos"
          + ("" if mal == 0 else "   <-- LA LEY NO CRUZA"))

print(SEP)
print(SEP)
if mal_train == 0 and total_mal == 0:
    print("RESULTADO: LEY CERRADA.  El signo en k=2 es la paridad de arriba, ajustada en m<=6 y")
    print("           validada sin un fallo en m=7,8,9, sin usar el senuelo.")
else:
    print("RESULTADO: NO.  La combinacion ajusta el entrenamiento pero no cruza a m=7,8,9:")
    print("           era un ajuste.  El signo sigue abierto.")
