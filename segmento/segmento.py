#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""segmento.py -- the conductor of a ring of character values, from the segment alone.

A single file, Python 3 standard library only.  No Sage, no numpy, no install.

WHAT IT COMPUTES
    Fix a torsion element of order q in SU(n), Sp(2m) or GL(n) and let A be the ring generated
    by the values of all characters at it.  This tool answers, for a prime p:

        familias(n, q)      which primes can divide [O:A] at all              (Thm. support)
        perfil(exps, M, p)  the layer dimensions W_0, W_1, ... of A at p      (Thm. local)
        indice(exps, M, p)  v_p([O:A]), as the sum of ALL the layer defects;
                            truncado(W, dim) says when that is a lower bound  (Thm. local (d))
        puente(n, q, p)     the same for the non-real ring of GL(n), at half
                            the work, via B = R (+) theta R                   (Thm. bridge)
        ceros(f)            how many components vanish on the diagonal,
                            f PRIME (it declines on composite f)     (note II, Cor. how many)
        fila(m)             whether row m is occupied, by the Jacobi symbol
                            (2|n) with n = 2m+1             (note II, Prop. rows, Thm. A)

WHY IT IS USEFUL OUTSIDE THIS PAPER
    The object is an order in a cyclotomic field and the question is its conductor.  That is the
    everyday computation in the algorithmics of complex multiplication -- endomorphism rings of
    abelian varieties, navigation of isogeny graphs -- where the conductor of a CM order is what
    one actually has to compute.  Two of the routines below replace a computation by a formula:

      * `familias` decides whether ANY prime divides the index, from divisibility alone: no
        lattice, no matrix.  It is a filter one can run over millions of pairs (n, q).
      * `puente` computes the index of the NON-REAL order from the real one.  It halves the
        degree: 2 v_p([O_+ : R]) + P d instead of closing a lattice in degree phi(q).

    Both are proved in the note (the first by the support theorem); `perfil` is the
    honest slow path and is used here as the control.

USAGE
    python segmento.py autotest              reproduce the worked examples of the notes
    python segmento.py perfil  n q p         the layer profile of GL(n) at p
    python segmento.py indice  n q p         v_p of the index of GL(n) (">=" when truncated)
    python segmento.py puente  n q p         the bridge shortcut, with its hypotheses checked
    python segmento.py familias n q          the primes that can divide the index
    python segmento.py ceros   f             vanishing components on the diagonal
    python segmento.py fila    m             can row m be occupied

    As a library:  from segmento import perfil, indice, puente, familias, ceros, fila

Authors: Carles Marin, Claude (AI assistant).  Public domain (CC0).
"""
import sys
from math import gcd

# --------------------------------------------------------------------------- aritmetica basica


def factoriza(n):
    """{primo: exponente}."""
    f, d = {}, 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def euler_phi(n):
    r = n
    for p in factoriza(n):
        r -= r // p
    return r


def valuacion(n, p):
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def orden_mult(a, m):
    """orden de a en (Z/m)^x."""
    a %= m
    if gcd(a, m) != 1:
        raise ValueError("a y m no son coprimos")
    k, x = 1, a
    while x != 1:
        x = (x * a) % m
        k += 1
    return k


def legendre2(n):
    """(2|n) para n impar, por la ley suplementaria: +1 si n = +-1 mod 8."""
    return 1 if n % 8 in (1, 7) else -1


# --------------------------------------------------------- polinomios sobre F_p, y ciclotomicos

def pol_recorta(a):
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def pol_mul(a, b, p):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] = (out[i + j] + x * y) % p
    return pol_recorta(out)


def pol_divmod(a, b, p):
    a = a[:]
    inv = pow(b[-1], p - 2, p)
    q = [0] * max(1, len(a) - len(b) + 1)
    while len(a) >= len(b) and any(a):
        d = len(a) - len(b)
        c = (a[-1] * inv) % p
        q[d] = c
        for i, y in enumerate(b):
            a[i + d] = (a[i + d] - c * y) % p
        pol_recorta(a)
        if len(a) < len(b):
            break
    return pol_recorta(q), pol_recorta(a)


def pol_mod(a, m, p):
    return pol_divmod(a, m, p)[1]


def ciclotomico(n):
    """Phi_n sobre Z, por division sucesiva de x^n - 1.  Coeficientes enteros exactos."""
    num = [-1] + [0] * (n - 1) + [1]
    for d in range(1, n):
        if n % d == 0:
            num = _div_exacta(num, ciclotomico(d))
    return num


_CICLO = {}


def _ciclo(n):
    if n not in _CICLO:
        _CICLO[n] = ciclotomico(n) if n > 1 else [-1, 1]
    return _CICLO[n]


def _div_exacta(a, b):
    """division exacta de polinomios sobre Z (b monico)."""
    a = a[:]
    q = [0] * (len(a) - len(b) + 1)
    for d in range(len(a) - len(b), -1, -1):
        c = a[d + len(b) - 1]
        q[d] = c
        if c:
            for i, y in enumerate(b):
                a[i + d] -= c * y
    return q


# ------------------------------------------------------------------- algebra lineal modulo p

def rango(filas, p):
    """rango de una matriz sobre F_p; no destruye la entrada."""
    M = [f[:] for f in filas]
    if not M:
        return 0
    nc = len(M[0])
    r = 0
    for c in range(nc):
        piv = None
        for i in range(r, len(M)):
            if M[i][c] % p:
                piv = i
                break
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], p - 2, p)
        M[r] = [(x * inv) % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] % p:
                f = M[i][c]
                M[i] = [(M[i][j] - f * M[r][j]) % p for j in range(nc)]
        r += 1
        if r == len(M):
            break
    return r


def _escalona(filas, p):
    """forma escalonada reducida, filas no nulas."""
    M = [f[:] for f in filas]
    nc = len(M[0]) if M else 0
    r = 0
    for c in range(nc):
        piv = None
        for i in range(r, len(M)):
            if M[i][c] % p:
                piv = i
                break
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], p - 2, p)
        M[r] = [(x * inv) % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] % p:
                f = M[i][c]
                M[i] = [(M[i][j] - f * M[r][j]) % p for j in range(nc)]
        r += 1
    return [f for f in M[:r]]


# --------------------------------------------------------------------------------- el segmento

def exps_SU(n):
    """exponentes del elemento principal de SU(n): n terminos de una clase de paridad mod 2q.
    Viven modulo 2q, no modulo q: usense con segmento_SU(n, q), que devuelve el par correcto."""
    return [n + 1 - 2 * j for j in range(1, n + 1)]


def segmento_SU(n, q):
    """(exponentes, M) del elemento principal g_q de SU(n), con la convencion de la nota:
    n par -> {n-1, n-3, ..., 1-n} modulo M = 2q (g_q tiene orden 2q);
    n impar -> los mismos exponentes divididos por 2, modulo M = q.
    Pasar exps_SU(n) con M = q para n par da OTRO elemento (el de exponentes impares mod q), y su
    anillo puede no ser el de g_q: SU(6) en q = 36 da v_3 = 1 asi y v_3 = 0 con g_q."""
    if n % 2 == 0:
        return exps_SU(n), 2 * q
    return [(n + 1 - 2 * j) // 2 for j in range(1, n + 1)], q


def exps_GL(n):
    return list(range(n))


def exps_Sp(m):
    out = []
    for j in range(1, m + 1):
        out += [j, -j]
    return out


def modulo_de(tipo, q, n=None):
    """el modulo M sobre el que viven los exponentes.  Para SU depende de la paridad de n (ver
    segmento_SU); hasta el 19-sep-2026 devolvia q siempre, lo que para n par es otro elemento."""
    if tipo in ("SU", "su"):
        if n is None:
            raise ValueError("SU: M depends on the parity of n; use segmento_SU(n, q)")
        return segmento_SU(n, q)[1]
    return q


# ------------------------------------------------------------------- las capas (el camino lento)

def perfil(exps, M, p, K=None):
    """dimensiones W_0, W_1, ... de la imagen del orden en O/(p, rad^K).

    Es el instrumento del paper: la imagen vive en F_p[y]/(Phi_r^K) con r la parte prima con p
    de M, y K no puede pasar de E = phi(p^v) porque mas alla la reduccion ya no es fiel.
    """
    v = valuacion(M, p)
    if v == 0:
        raise ValueError("p no divide a M: no hay nada que medir")
    r = M // p ** v
    E = euler_phi(p ** v)
    K = E if K is None else min(K, E)
    phi = _ciclo(r)
    phi_p = [c % p for c in phi]
    mod = [1]
    for _ in range(K):
        mod = pol_mul(mod, phi_p, p)
    D = len(mod) - 1

    def red(a):
        return pol_mod(a, mod, p)

    # polinomio caracteristico por producto en arbol
    facs = []
    for e in exps:
        yk = [0] * ((e % M) + 1)
        yk[e % M] = 1
        facs.append([red([(-x) % p for x in yk]), [1]])   # (T - y^e) -> [-y^e, 1]
    while len(facs) > 1:
        nue = []
        for i in range(0, len(facs) - 1, 2):
            a, b = facs[i], facs[i + 1]
            out = [[0]] * (len(a) + len(b) - 1)
            out = [[0] for _ in range(len(a) + len(b) - 1)]
            for ia, ca in enumerate(a):
                if ca == [0]:
                    continue
                for ib, cb in enumerate(b):
                    if cb == [0]:
                        continue
                    t = pol_mul(ca, cb, p)
                    s = out[ia + ib]
                    ln = max(len(s), len(t))
                    s = s + [0] * (ln - len(s))
                    t = t + [0] * (ln - len(t))
                    out[ia + ib] = pol_recorta([(s[i] + t[i]) % p for i in range(ln)])
            facs.append(out) if False else None
            nue.append([red(c) for c in out])
        if len(facs) % 2:
            nue.append(facs[-1])
        facs = nue
    co = facs[0]
    nn = len(exps)
    gens = []
    for k in range(1, nn + 1):
        c = co[nn - k]
        g = [((-1) ** k * x) % p for x in c]
        g = pol_recorta(g)
        if g != [0]:
            gens.append(g)

    def vec(a):
        return (a + [0] * (D - len(a)))[:D]

    base = _escalona([vec([1])] + [vec(g) for g in gens], p)
    while True:
        prods = []
        for b in base:
            eb = pol_recorta(b[:])
            for g in gens:
                prods.append(vec(red(pol_mul(eb, g, p))))
        nueva = _escalona(base + prods, p)
        if len(nueva) == len(base):
            break
        base = nueva

    dims = []
    for k in range(1, K + 1):
        mk = [1]
        for _ in range(k):
            mk = pol_mul(mk, phi_p, p)
        Dk = len(mk) - 1
        filas = [(pol_mod(pol_recorta(b[:]), mk, p) + [0] * Dk)[:Dk] for b in base]
        dims.append(rango(filas, p))
    return [dims[0]] + [dims[i] - dims[i - 1] for i in range(1, len(dims))]


def indice(exps, M, p, dim=None):
    """v_p([O:A]) = suma de (dim - W_i) sobre las capas defectuosas."""
    v = valuacion(M, p)
    r = M // p ** v
    if r <= 2:
        # GUARDA: qprima<=2
        # Fuera del teorema local: con q' <= 2 el algebra residual es F_p, d = phi(q')/2 no es la
        # dimension y la matriz de momentos no tiene columnas.  Hasta el 19-sep-2026 esta rutina
        # devolvia v_p = 0 aqui, y SU(4) en q = 8, p = 2 tiene indice 2.
        raise ValueError("q' = %d <= 2: outside the local theorem, the index is not computed here" % r)
    W = perfil(exps, M, p)
    if dim is None:
        dim = euler_phi(r)
        if _es_simetrico(exps, M):
            dim //= 2
    # Se suman TODAS las capas.  Hasta el 19-sep-2026 se paraba en la primera capa llena, que es
    # suponer "una capa llena => las siguientes tambien" -- el paso que la nota declara NO
    # demostrado -- y en SU(15), q=45, p=3 el perfil [1,0,1,2,1,2] vuelve a bajar tras llenarse:
    # la longitud es 5, no 4.
    total = sum(max(0, dim - w) for w in W)
    return total, W, dim


def truncado(W, dim):
    """True si el perfil se corta en profundidad K = E ANTES de llenarse: entonces la suma de los
    defectos es una COTA INFERIOR de v_p, no su valor.  Ejemplo: GL(15), q=45, p=3 da
    [1,0,1,3,2,3] con dimension 4 y la ultima capa no esta llena."""
    return bool(W) and W[-1] < dim


def _es_simetrico(exps, M):
    s = sorted(e % M for e in exps)
    t = sorted((-e) % M for e in exps)
    return s == t


# ------------------------------------------------------------------------------- las formulas

def familias(n, q):
    """primos que PUEDEN dividir el indice: aquellos con r = q/p^{v_p(q)} en una familia.

    Es la forma util del criterio de soporte: sin calcular nada, descarta.
    """
    out = []
    for p in sorted(factoriza(q)):
        r = q // p ** valuacion(q, p)
        if any((n - 1) % r == 0 or n % r == 0 or (n + 1) % r == 0 for _ in (0,)):
            out.append(p)
    return out


def puente(n, q, p):
    """el atajo del puente: v_p([O:B]) = 2 v_p([O_+ : R]) + P d, con P = p^{v_p(N)}.

    Devuelve (v_p, P, d, hipotesis) donde `hipotesis` dice si el atajo aplica.
    """
    v = valuacion(q, p)
    r = q // p ** v
    hip = []
    if p == 2:
        hip.append("p = 2: the shortcut does not exist, the proof divides by 2")
    if r % 2 == 0:
        hip.append("r even: outside the family treated")
    if r < 5:
        hip.append("r < 5")
    if n % r == 0:
        N = n // r
    elif (n - 1) % r == 0:
        N = (n - 1) // r
    else:
        hip.append("n is neither Nr nor Nr+1: outside the family")
        N = None
    if N is None or not (2 <= n < q):
        return None, None, None, hip
    P = p ** valuacion(N, p)
    if P >= p ** v:
        hip.append("P >= p^v: the segment closes, the extension degenerates")
    d = euler_phi(r) // 2
    # el lado real: el segmento recentrado, simetrico
    inv2 = pow(2, -1, q)
    a = (N * r * inv2) % q
    S = list(range(1, N * r)) if n % r == 0 else list(range(0, N * r + 1))
    cent = [(e - a) % q for e in S]
    vR, WR, dR = indice(cent, q, p, dim=d)
    if truncado(WR, d):
        # GUARDA: truncado
        hip.append("the real profile is cut at depth E before it fills: the value is a LOWER BOUND")
    # El atajo necesita R = A.  Hasta el 18-sep-2026 eso estaba probado solo bajo u > d/2 y aqui
    # habia una guarda; prop:tau lo demuestra ahora SIN condiciones (dos identidades en
    # Z[t,t^-1] con denominador unidad), asi que la guarda se retiro el 19-sep: declinar un
    # numero justificado es tan falso como devolver uno que no lo esta.
    return 2 * vR + P * d, P, d, hip


def _es_primo(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def ceros(f):
    """cuantas componentes se anulan en la diagonal, y cuales caracteres.

    0 si -1 pertenece a <2>, y d/t si no, con t el orden de 2 modulo f.

    EXIGE f PRIMO, y no es un capricho.  La formula sale de que B_{1,chi} no se anula, y eso lo
    da la primitividad; para f compuesto un caracter impar IMPRIMITIVO puede tener B_{1,chi} = 0
    por sus factores de Euler, y entonces la componente se anula sin que chi(2) = 1.  El fallo
    mas pequeno es f = 21: la formula predice una componente nula y hay dos.  Enumerando todo
    f <= 121 impar, acierta en los 28 primos y falla en 12 de los 31 compuestos.

    Devolver d/t para f compuesto seria devolver un numero equivocado con cara de correcto, que
    es peor que no devolver nada.
    """
    if f % 2 == 0 or f < 5:
        raise ValueError("f debe ser impar y al menos 5")
    if not _es_primo(f):
        # GUARDA: f-compuesto
        raise ValueError(
            "f = %d es compuesto: la formula d/t solo vale para f primo.  Para f compuesto un "
            "caracter impar imprimitivo puede anular B_{1,chi} y la cuenta falla (el primer caso "
            "es f = 21, donde predice 1 y hay 2).  Esta rutina declina en vez de devolver un "
            "numero que no puede justificar." % f)
    d = euler_phi(f) // 2
    t = orden_mult(2, f)
    pot = {pow(2, j, f) for j in range(t)}
    if (f - 1) in pot:
        return 0, t, d
    return d // t, t, d


def fila(m):
    """puede la fila m tener ceros?  Solo si el SIMBOLO DE JACOBI (2|n) vale 1, con n = 2m+1.

    NO es "2 es resto cuadratico modulo n", que era como estaba escrito y es otra afirmacion,
    mas fuerte: para n compuesto el simbolo de Jacobi puede valer +1 sin que 2 sea un cuadrado.
    En n = 15 se tiene (2|15) = +1 y los cuadrados modulo 15 son 0,1,4,6,9,10, sin el 2.  Las dos
    lecturas coinciden solo para n primo.  Lo que el criterio dice es (2|n) = 1, equivalentemente
    n = +-1 modulo 8, y es lo que legendre2() calcula.

    Las dos mitades estan demostradas:
      (2|n) = -1  =>  la fila esta VACIA: la suma pesada tiene la paridad del numero triangular
                      T_m = (n^2-1)/8, y (2|n) = (-1)^{T_m}  (nota II, Prop. rows).
      (2|n) = +1  =>  la fila esta OCUPADA: existe una solucion completamente multiplicativa
                      (nota II, Thm. every admissible length, y Cor. rows).
    Devuelve (ocupada, n, T_m, estatuto).
    """
    n = 2 * m + 1
    T = m * (m + 1) // 2
    if legendre2(n) == -1:
        return False, n, T, "proved: the row is empty (note II, Prop. rows)"
    return True, n, T, "proved: the row is occupied (note II, Thm. every admissible length)"


# ------------------------------------------------------------------------------- autoprueba

CASOS_PERFIL = [
    # (exponentes, M, p, perfil esperado)  -- cifras que aparecen en la nota
    (exps_GL(23), 207, 3, [1, 11, 21, 22]),
    (exps_GL(24), 207, 3, [1, 11, 21, 22]),
    (exps_GL(31), 279, 3, [1, 15, 29, 30]),
    (exps_GL(11), 99, 3, [1, 6, 10, 10]),
    (exps_GL(15), 45, 3, [1, 0, 1, 3, 2, 3]),
]

CASOS_INDICE = [
    (exps_GL(23), 207, 3, 33),
    (exps_GL(11), 99, 3, 13),
    (exps_GL(7), 63, 3, 7),
    (exps_GL(5), 45, 3, 4),
]

CASOS_CEROS = [(7, 1), (23, 1), (31, 3), (47, 1), (11, 0), (13, 0), (29, 0)]
# Los compuestos donde la formula FALLA: la rutina tiene que declinar, no devolver d/t.  Sin
# estos casos el autotest solo probaba la mitad que no puede fallar -- todos los f de arriba son
# primos -- y una guarda que nadie ejercita es una guarda que nadie sabe si existe.
CASOS_CEROS_DECLINA = [21, 33, 39, 55, 57, 63, 93, 95, 99, 105, 111, 117]
CASOS_FILA = [(3, True), (4, True), (5, False), (6, False), (7, True), (15, True), (17, False)]


def autotest():
    fallos = 0
    print("layer profile")
    for exps, M, p, esp in CASOS_PERFIL:
        got = perfil(exps, M, p)[:len(esp)]
        ok = got == esp
        fallos += 0 if ok else 1
        print("   M=%-5d p=%-3d  %-26s expected %-26s %s" % (
            M, p, got, esp, "ok" if ok else "*** FAILS ***"))
    print("index v_p")
    for exps, M, p, esp in CASOS_INDICE:
        got = indice(exps, M, p)[0]
        ok = got == esp
        fallos += 0 if ok else 1
        print("   M=%-5d p=%-3d  v_p=%-4d expected %-4d %s" % (
            M, p, got, esp, "ok" if ok else "*** FAILS ***"))
    print("the bridge shortcut against the slow path")
    for n, q, p in [(23, 207, 3), (11, 99, 3), (7, 63, 3), (5, 45, 3), (15, 45, 3)]:
        rapido, P, d, hip = puente(n, q, p)
        dimB = euler_phi(q // p ** valuacion(q, p))
        lento, Wl, _ = indice(exps_GL(n), q, p, dim=dimB)
        cota = truncado(Wl, dimB)
        ok = (rapido == lento)
        fallos += 0 if ok else 1
        print("   n=%-4d q=%-5d p=%-3d  short=%-5s slow=%s%-5d P=%-3s %s%s" % (
            n, q, p, rapido, ">=" if cota else "", lento, P, "ok" if ok else "*** FAILS ***",
            "  (both are lower bounds: the profile is cut before it fills)" if cota else ""))
    print("vanishing components on the diagonal")
    for f, esp in CASOS_CEROS:
        got = ceros(f)[0]
        ok = got == esp
        fallos += 0 if ok else 1
        print("   f=%-5d  %d expected %d %s" % (f, got, esp, "ok" if ok else "*** FAILS ***"))
    print("composite moduli: the routine must decline, not answer")
    for f in CASOS_CEROS_DECLINA:
        try:
            ceros(f)
            declina = False
        except ValueError:
            declina = True
        fallos += 0 if declina else 1
        print("   f=%-5d  %s" % (f, "declines, as it must" if declina
                                 else "*** FAILS: answered instead of declining ***"))
    print("rows that can be occupied")
    for m, esp in CASOS_FILA:
        got = fila(m)[0]
        ok = got == esp
        fallos += 0 if ok else 1
        print("   m=%-4d  %-6s expected %-6s %s" % (m, got, esp, "ok" if ok else "*** FAILS ***"))
    print()
    print("SELF-TEST: %s" % ("all correct" if fallos == 0 else "%d FAILURES" % fallos))
    return fallos


# ------------------------------------------------------------------------------------- la CLI

def main(argv):
    if len(argv) < 2 or argv[1] in ("-h", "--help"):
        print(__doc__)
        return 0
    cmd = argv[1]
    if cmd == "autotest":
        return 1 if autotest() else 0
    if cmd == "familias":
        n, q = int(argv[2]), int(argv[3])
        ps = familias(n, q)
        print("primes that can divide the index of (n=%d, q=%d): %s" % (
            n, q, ps if ps else "none -- the order is maximal (the support theorem)"))
        return 0
    if cmd == "ceros":
        f = int(argv[2])
        c, t, d = ceros(f)
        print("f=%d: order of 2 = %d, d = %d, vanishing components = %d" % (f, t, d, c))
        return 0
    if cmd == "fila":
        m = int(argv[2])
        ok, n, T, est = fila(m)
        print("m=%d: n = %d, T_m = %d, (2|n) = %s -> the row %s" % (
            m, n, T, legendre2(n), "is occupied" if ok else "is empty"))
        print("   status: %s" % est)
        return 0
    n, q, p = int(argv[2]), int(argv[3]), int(argv[4])
    if cmd == "perfil":
        print(perfil(exps_GL(n), q, p))
        return 0
    if cmd == "indice":
        v, W, dim = indice(exps_GL(n), q, p)
        print("v_p %s %d   (layers %s, full dimension %d)%s" % (
            ">=" if truncado(W, dim) else "=", v, W, dim,
            "   LOWER BOUND: the profile is cut at depth E before it fills" if truncado(W, dim) else ""))
        return 0
    if cmd == "puente":
        val, P, d, hip = puente(n, q, p)
        if hip:
            print("hypotheses not met: " + "; ".join(hip))
        if val is not None:
            print("v_p([O:B]) = 2 v_p([O_+:R]) + P d = %d   with P = %d, d = %d" % (val, P, d))
        return 0
    print("unknown command: %s" % cmd)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv))
