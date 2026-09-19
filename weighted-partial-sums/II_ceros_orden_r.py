# -*- coding: utf-8 -*-
r"""II_ceros_orden_r.py -- los ceros (q^j-1)/2 para caracteres impares de CUALQUIER orden con chi(2)=1.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

chi: caracter de Dirichlet mod q primo, extendido completamente multiplicativo con chi(q) = u
(u raiz de la unidad cualquiera).  Afirmacion: si chi es impar y chi(2) = 1, entonces
S((q^j - 1)/2) = sum_{n <= (q^j-1)/2} n chi(n) = 0 para todo j y todo u.
(Prueba: (D) vale igual; T((q-1)/2) = 0 por el teorema de la diagonal; C((q-1)/2) = U =
(chi(2)^{-1} - 2) B_{1,chi} = -B_{1,chi} y T_q = q B_{1,chi}, asi que T_q + q C((q-1)/2) = 0.)
Calculo EXACTO: chi(n) = zeta_r^{ind(n)}, sumas en Z[zeta_r] como vectores de enteros por exponente
(se reduce al final: S = 0 en Z[zeta_r] sii el vector es combinacion de las relaciones de Phi_r;
se comprueba reduciendo el polinomio mod Phi_r con enteros).
Para cada primo q <= QMAX y cada caracter impar (orden r | q-1, par) con chi(2) = 1: j = 1, 2 (y 3 si
cabe bajo JMAX), u en {1, zeta_r} (y u = -1).  Control: los impares con chi(2) != 1 NO se anulan en j=1.
Cuenta ademas cuantos caracteres impares tienen chi(2) = 1 frente a d/t (Cor. de la diagonal).
Uso: python II_ceros_orden_r.py [QMAX] [JMAX_N]
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 130
NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 3_000_000


def es_primo(n):
    return n > 1 and all(n % d for d in range(2, int(n ** .5) + 1))


def raiz_primitiva(q):
    fs = [p for p in range(2, q) if (q - 1) % p == 0 and es_primo(p)]
    for g in range(2, q):
        if all(pow(g, (q - 1) // p, q) != 1 for p in fs):
            return g


def cicl(r):
    """coeficientes enteros de Phi_r (grado phi(r)), por division de x^r - 1."""
    import functools
    def mul(a, b):
        c = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                c[i + j] += x * y
        return c
    def div(a, b):
        a = a[:]; q = [0] * (len(a) - len(b) + 1)
        for i in range(len(a) - len(b), -1, -1):
            c = a[i + len(b) - 1] // b[-1]; q[i] = c
            for j, y in enumerate(b):
                a[i + j] -= c * y
        return q
    P = [-1] + [0] * (r - 1) + [1]
    for d in range(1, r):
        if r % d == 0:
            P = div(P, cicl(d))
    return P


def es_cero(vec, r):
    """vec[k] = coeficiente de zeta^k (k < r); cero en Z[zeta_r] sii Phi_r divide el polinomio."""
    Phi = cicl(r)
    a = vec[:]
    for i in range(len(a) - 1, len(Phi) - 2, -1):
        c = a[i]
        if c:
            for j, y in enumerate(Phi):
                a[i - (len(Phi) - 1) + j] -= c * y
    return all(x == 0 for x in a)


def S_exacto(q, g, e, r, uexp, m):
    """sum_{n<=m} n chi(n), chi(g) = zeta_r^e, chi(q) = zeta_r^uexp; vector por exponente."""
    ind = [0] * q
    x = 1
    for k in range(q - 1):
        ind[x] = k; x = x * g % q
    vec = [0] * r
    for n in range(1, m + 1):
        v = 0; y = n
        while y % q == 0:
            y //= q; v += 1
        k = (e * ind[y % q] + uexp * v) % r
        vec[k] += n
    return vec


tot = fallos = ctrl = ctrl_mal = 0
cuenta_ok = 0; cuenta_tot = 0
for q in range(5, QMAX + 1):
    if not es_primo(q):
        continue
    g = raiz_primitiva(q)
    t = 1; x = 2
    while x != 1:
        x = x * 2 % q; t += 1
    d = (q - 1) // 2
    impares_2 = 0
    for e in range(1, q - 1):
        # caracter chi(g) = zeta_{q-1}^e; orden r = (q-1)/gcd
        from math import gcd
        r = (q - 1) // gcd(e, q - 1)
        ee = e // gcd(e, q - 1)
        # impar: chi(-1) = chi(g)^{(q-1)/2} = -1  <=>  ee*(q-1)/2 ... en orden r: exponente ee*(r/2) si r par
        # chi(-1) = chi(g)^{(q-1)/2} = zeta_r^{ee (q-1)/2}; impar sii ese exponente es r/2 mod r
        if r % 2 or (ee * ((q - 1) // 2)) % r != r // 2:
            continue
        # chi(2): 2 = g^ind2
        ind2 = next(k for k in range(q - 1) if pow(g, k, q) == 2)
        chi2_1 = (ee * ind2) % r == 0
        if chi2_1:
            impares_2 += 1
            for uexp in sorted({0, 1 % r, (r // 2) % r}):
                for j in (1, 2, 3):
                    m = (q ** j - 1) // 2
                    if m > NMAX:
                        break
                    tot += 1
                    if not es_cero(S_exacto(q, g, ee, r, uexp, m), r):
                        fallos += 1
                        print("FALLA q=%d r=%d e=%d u=zeta^%d j=%d" % (q, r, ee, uexp, j))
        else:
            ctrl += 1
            if es_cero(S_exacto(q, g, ee, r, 0, (q - 1) // 2), r):
                ctrl_mal += 1
    esperado = 0 if pow(2, t // 2, q) == q - 1 and t % 2 == 0 else d // t
    # -1 en <2> sii t par y 2^{t/2} = -1
    cuenta_tot += 1; cuenta_ok += (impares_2 == esperado)
print("primos q <= %d; comprobaciones (chi impar, chi(2)=1, u, j): %d, fallos: %d" % (QMAX, tot, fallos))
print("control: impares con chi(2) != 1 que se anulan en j=1: %d de %d" % (ctrl_mal, ctrl))
print("numero de impares con chi(2)=1 igual a d/t (o 0): %d de %d primos" % (cuenta_ok, cuenta_tot))
