# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
#
# QUE: nucleo de calculo (python puro, enteros exactos) para ordenes del tipo
#        A = Z[ e_k(eta_a : a en X+) ]  dentro de  O_N = Z[xi + xi^-1],  eta_a = xi^a + xi^-a,  xi = zeta_N,
#      para una lista cualquiera X+ de exponentes positivos (segmento lleno {1..m} del articulo, o el segmento
#      impar {1,3,...,2m-1} de los tipos A_{2m-1} y B_m).  Calcula:
#        - el reticulo de A en la base theta^k (theta = xi + xi^-1) por cierre multiplicativo con HNF (modular
#          en cuanto tiene rango lleno), el indice [O:A] y su factorizacion;
#        - la estabilidad de Galois de A (sigma_u(e_k) en A, via Chebyshev);
#        - en cada primo p | N, con N = p^v N': las capas W_i del articulo (filtracion por rad^i, rad = z O_p con
#          z = psi_{N'}(theta), uniformizante en todos los primos sobre p) y el exponente e_p del conductor
#          (primer i con rad^i contenido en A_p);
#        - la prediccion del primer momento: S(c) = suma de x en X = +-X+, x = c (mod N'), y
#          rango_{F_p} [S(u^{-1} c)]_{u, c};  la estabilidad del multiconjunto reducido mod N'.
# POR QUE: falsar si el enlace "primer momento = Stickelberger => numeros de clase deciden el conductor" es ley
#      mas alla del segmento lleno {1..m} (tipo C).  Es el mismo algoritmo que la nota (reticulo + capas), escrito
#      de nuevo en python puro para no depender de Sage.
# FECHA: 28-sep-2026
from math import gcd


def phi(n):
    r, x, p = n, n, 2
    while p * p <= x:
        if x % p == 0:
            while x % p == 0:
                x //= p
            r -= r // p
        p += 1
    if x > 1:
        r -= r // x
    return r


def factor(n):
    f, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def vp(x, p):
    v = 0
    while x and x % p == 0:
        x //= p
        v += 1
    return v


# ---------------------------------------------------------------- polinomios enteros
def pmul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] += x * y
    return r


def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]


def cyclo(n):
    """Phi_n(x) con coeficientes enteros (lista, grado creciente)."""
    num = [-1] + [0] * (n - 1) + [1]  # x^n - 1
    for dd in range(1, n):
        if n % dd == 0:
            num = pdiv_exact(num, cyclo(dd))
    return num


_cyc = {}


def pdiv_exact(a, b):
    a = a[:]
    q = [0] * (len(a) - len(b) + 1)
    for i in range(len(a) - len(b), -1, -1):
        c = a[i + len(b) - 1] // b[-1]
        q[i] = c
        for j, y in enumerate(b):
            a[i + j] -= c * y
    assert all(x == 0 for x in a), "division no exacta"
    return q


def cheb_polys(n):
    """C_0..C_n como polinomios en T: C_u(t + 1/t) = t^u + t^-u (C_0 = 2)."""
    C = [[2], [0, 1]]
    for u in range(1, n):
        C.append(padd(pmul([0, 1], C[u]), [-x for x in C[u - 1]]))
    return C[: n + 1]


def psi(n):
    """Polinomio minimo (monico, entero) de 2cos(2 pi / n); psi_1 = T-2, psi_2 = T+2."""
    if n == 1:
        return [-2, 1]
    if n == 2:
        return [2, 1]
    f = cyclo(n)
    d = (len(f) - 1) // 2
    C = cheb_polys(d)
    r = [f[d]]
    for k in range(1, d + 1):
        r = padd(r, [f[d + k] * x for x in C[k]])
    while len(r) > 1 and r[-1] == 0:
        r.pop()
    assert r[-1] == 1 and len(r) - 1 == d
    return r


class Anillo:
    """O_N = Z[T]/psi_N(T); elementos = listas de d enteros."""

    def __init__(self, N):
        self.N = N
        self.f = psi(N)
        self.d = len(self.f) - 1

    def red(self, a, mod=None):
        a = list(a)
        d, f = self.d, self.f
        for i in range(len(a) - 1, d - 1, -1):
            c = a[i]
            if c:
                for j in range(d + 1):
                    a[i - d + j] -= c * f[j]
        a = (a + [0] * d)[:d]
        if mod:
            a = [x % mod for x in a]
        return a

    def mul(self, a, b, mod=None):
        return self.red(pmul(a, b), mod)

    def uno(self):
        return [1] + [0] * (self.d - 1)

    def evalpoly(self, P, x, mod=None):
        """P(x) con P entero en T y x en O (Horner)."""
        r = [0] * self.d
        for c in reversed(P):
            r = self.mul(r, x, mod) if any(r) else [0] * self.d
            r[0] += c
            if mod:
                r = [y % mod for y in r]
        return r


def generadores(O, exps):
    """e_1..e_m de los eta_a, a en exps (en la base theta^k)."""
    C = cheb_polys(max(exps))
    etas = [O.red(C[a]) for a in exps]
    coef = [O.uno()]  # prod (Y - eta), coeficientes de Y^0.. en O
    for eta in etas:
        new = [[0] * O.d for _ in range(len(coef) + 1)]
        for i, c in enumerate(coef):
            new[i + 1] = [x + y for x, y in zip(new[i + 1], c)]
            ec = O.mul(eta, c)
            new[i] = [x - y for x, y in zip(new[i], ec)]
        coef = new
    m = len(exps)
    # coef[m-k] = (-1)^k e_k
    return [[(-1) ** k * x for x in coef[m - k]] for k in range(1, m + 1)], etas


# ---------------------------------------------------------------- reticulos (HNF por insercion)
def egcd(a, b):
    x0, y0, x1, y1 = 1, 0, 0, 1
    while b:
        q, a, b = a // b, b, a % b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return a, x0, y0


class Reticulo:
    def __init__(self, d, D=None):
        self.d = d
        self.piv = {}
        self.D = D  # si no es None: D Z^d esta contenido en el reticulo

    def _norm(self, v):
        if self.D:
            return [x % self.D for x in v]
        return v

    def insertar(self, v):
        v = self._norm(list(v))
        cambio = False
        for col in range(self.d):
            b = v[col]
            if self.D:
                b %= self.D
                v[col] = b
            if b == 0:
                continue
            if col not in self.piv:
                if self.D:
                    g, x, _ = egcd(b, self.D)
                    self.piv[col] = self._norm([x * y for y in v])
                    # el resto (D/g) v tiene 0 en col
                    v = self._norm([(self.D // g) * y for y in v])
                    v[col] = 0
                    cambio = True
                    continue
                if b < 0:
                    v = [-y for y in v]
                self.piv[col] = v
                return True
            P = self.piv[col]
            a = P[col]
            if b % a == 0:
                q = b // a
                v = self._norm([y - q * z for y, z in zip(v, P)])
                continue
            g, x, y = egcd(a, b)
            if g < 0:
                g, x, y = -g, -x, -y
            new = self._norm([x * s + y * t for s, t in zip(P, v)])
            rem = self._norm([(a // g) * t - (b // g) * s for s, t in zip(P, v)])
            self.piv[col] = new
            v = rem
            cambio = True
        if self.D and cambio:
            self._actualizar_D()
        return cambio

    def _actualizar_D(self):
        if len(self.piv) == self.d:
            det = 1
            for c in range(self.d):
                det *= self.piv[c][c]
            if det < self.D:
                self.D = det
                self._reducir_filas()

    def _reducir_filas(self):
        """filas pivote mod D, conservando la diagonal (puede valer D si las demas son 1)."""
        for c in self.piv:
            a = self.piv[c][c]
            fila = [x % self.D for x in self.piv[c]]
            fila[c] = a
            self.piv[c] = fila

    def rango(self):
        return len(self.piv)

    def indice(self):
        """[Z^d : L] si rango lleno (en modo modular faltan pivotes = D)."""
        det = 1
        for c in range(self.d):
            if c in self.piv:
                det *= self.piv[c][c] if not self.D else gcd(self.piv[c][c], self.D)
            else:
                if not self.D:
                    return None
                det *= self.D
        return det

    def contiene(self, v):
        v = self._norm(list(v))
        for col in range(self.d):
            b = v[col] % self.D if self.D else v[col]
            if b == 0:
                continue
            if col not in self.piv:
                return False
            a = self.piv[col][col]
            if b % a:
                return False
            q = b // a
            v = self._norm([y - q * z for y, z in zip(v, self.piv[col])])
        return True

    def modular(self):
        """Pasa a modo modular con D = det (requiere rango lleno)."""
        if self.D or len(self.piv) < self.d:
            return
        det = 1
        for c in range(self.d):
            det *= self.piv[c][c]
        self.D = det
        self._reducir_filas()


def cierre(O, gens):
    """Reticulo del subanillo Z[gens] de O (base theta^k)."""
    L = Reticulo(O.d)
    cola = []
    for v in [O.uno()] + gens:
        if L.insertar(v):
            cola.append(v)
    while cola:
        v = cola.pop()
        for g in gens:
            w = O.mul(g, v, L.D)
            if L.insertar(w):
                cola.append(w)
                if not L.D and L.rango() == O.d:
                    L.modular()
    if L.rango() == O.d:
        L.modular()
    return L


# ---------------------------------------------------------------- orden, Galois, capas
def unidades_mod_pm(n):
    """representantes de (Z/n)^x / +-1 (el menor de u, n-u)."""
    if n <= 2:
        return [1]
    return [u for u in range(1, n // 2 + 1) if gcd(u, n) == 1 and u <= n - u]


class Orden:
    def __init__(self, N, exps):
        self.N = N
        self.exps = list(exps)
        self.O = Anillo(N)
        self.d = self.O.d
        self.gens, self.etas = generadores(self.O, self.exps)
        self.L = cierre(self.O, self.gens)
        self.rango = self.L.rango()
        self.idx = self.L.indice() if self.rango == self.d else None

    def galois_estable(self):
        """sigma_u(e_k) en A para toda unidad u; theta -> C_u(theta)."""
        if self.idx is None:
            return None
        C = cheb_polys(self.N)
        for u in unidades_mod_pm(self.N):
            if u == 1:
                continue
            su = self.O.red(C[u])
            for g in self.gens:
                img = self.O.evalpoly(g, su, self.L.D)
                if not self.L.contiene(img):
                    return False
        return True

    def capas(self, p, extra=2, imax=None):
        """(W, e, a): capas W_0..W_{e+extra-1} y exponente e del conductor en p.  a = v_p[O:A]."""
        a = vp(self.idx, p)
        if a == 0:
            return [], 0, 0, []
        v = vp(self.N, p)
        Np = self.N // p ** v
        z = self.O.red(psi(Np))
        mod = p ** a
        O = self.O
        base_A = [self.L.piv[c] for c in range(self.d)]

        def ell(zi):
            """(v_p[O : A + zi O + p^a O], v_p[O : zi O + p^a O])."""
            J = Reticulo(self.d, mod)
            gensJ = []
            x = [y % mod for y in zi]
            t = O.uno()
            for k in range(self.d):
                gensJ.append(O.mul(x, t, mod))
                t = O.mul(t, [0, 1] + [0] * (self.d - 2), mod) if self.d > 1 else t
            for g in gensJ:
                J.insertar(g)
            nJ = vp(J.indice(), p)
            AJ = Reticulo(self.d, mod)
            for g in gensJ + base_A:
                AJ.insertar(g)
            return vp(AJ.indice(), p), nJ

        ls, ns = [], []
        zi = O.uno()
        i = 0
        e = None
        tope = imax if imax else a * phi(p ** v) * (2 if Np <= 2 else 1) + 1
        while True:
            l_i, n_i = ell(zi)
            ls.append(l_i)
            ns.append(n_i)
            if e is None and l_i == a:
                e = i
            if e is not None and i >= e + extra:
                break
            if i >= tope:
                break
            zi = O.mul(zi, z, mod)
            i += 1
        W = [(ns[k + 1] - ls[k + 1]) - (ns[k] - ls[k]) for k in range(len(ls) - 1)]
        dres = [ns[k + 1] - ns[k] for k in range(len(ns) - 1)]
        return W, e, a, dres


def capas_mod(N, exps, p, t=1, imax=None):
    """Capas W_0..W_{imax-1} de A en p calculadas SOLO modulo p^t (sin reticulo global): A + p^t O por cierre
    modular.  Validas para i < t*E (rad^{tE} = p^t O_p).  Devuelve (W, dres)."""
    O = Anillo(N)
    gens, _ = generadores(O, exps)
    mod = p ** t
    L = Reticulo(O.d, mod)
    cola = []
    for v in [O.uno()] + gens:
        if L.insertar(v):
            cola.append(v)
    while cola:
        v = cola.pop()
        for g in gens:
            w = O.mul(g, v, mod)
            if L.insertar(w):
                cola.append(w)
    v = vp(N, p)
    Np = N // p ** v
    z = O.red(psi(Np), mod)
    E = phi(p ** v)
    tope = imax if imax else t * E
    base_A = [L.piv[c] for c in sorted(L.piv)]
    theta = [0, 1] + [0] * (O.d - 2)
    ls, ns = [], []
    zi = O.uno()
    for i in range(tope + 1):
        J = Reticulo(O.d, mod)
        AJ = Reticulo(O.d, mod)
        tk = [y % mod for y in zi]
        gensJ = []
        for k in range(O.d):
            gensJ.append(tk)
            tk = O.mul(tk, theta, mod)
        for g in gensJ:
            J.insertar(g)
        for g in gensJ + base_A:
            AJ.insertar(g)
        ns.append(vp(J.indice(), p))
        ls.append(vp(AJ.indice(), p))
        zi = O.mul(zi, z, mod)
    W = [(ns[k + 1] - ls[k + 1]) - (ns[k] - ls[k]) for k in range(tope)]
    dres = [ns[k + 1] - ns[k] for k in range(tope)]
    return W, dres


def local_completo(N, exps, p, tmax=40):
    """Sin reticulo global: a = v_p[O:A] por cierre mod p^t con t creciente hasta que v_p[O : A + p^t O] se
    estabiliza (entonces p^t O esta en A_p, Nakayama), y capas W_0..W_{(a+1)E-1} con t = a+1; e = primer i desde el
    que todas las capas calculadas son llenas (valido: rad^{aE} = p^a O_p esta en A_p).  Devuelve (a, e, W)."""
    O = Anillo(N)
    gens, _ = generadores(O, exps)
    prev = None
    for t in range(1, tmax + 1):
        mod = p ** t
        L = Reticulo(O.d, mod)
        cola = []
        for v in [O.uno()] + gens:
            if L.insertar(v):
                cola.append(v)
        while cola:
            v = cola.pop()
            for g in gens:
                w = O.mul(g, v, mod)
                if L.insertar(w):
                    cola.append(w)
        at = vp(L.indice(), p)
        if prev is not None and at == prev:
            a = at
            break
        prev = at
    else:
        raise RuntimeError("no se estabiliza")
    W, dres = capas_mod(N, exps, p, a + 1)
    d = dres[0]
    e = len(W)
    while e > 0 and W[e - 1] == d:
        e -= 1
    return a, e, W


# ---------------------------------------------------------------- primer momento
def conjunto_signado(exps):
    return [x for a in exps for x in (a, -a)]


def S_de(exps, n):
    S = [0] * n
    for x in conjunto_signado(exps):
        S[x % n] += x
    return S


def cuentas(exps, n):
    c = [0] * n
    for x in conjunto_signado(exps):
        c[x % n] += 1
    return c


def estabilizador(exps, n):
    """u en (Z/n)^x/+-1 que fijan el multiconjunto reducido."""
    c = cuentas(exps, n)
    return [u for u in unidades_mod_pm(n) if all(c[(u * x) % n] == c[x] for x in range(n))]


def rango_mod(filas, p):
    M = [[x % p for x in f] for f in filas]
    r = 0
    ncols = len(M[0]) if M else 0
    for col in range(ncols):
        pv = next((i for i in range(r, len(M)) if M[i][col]), None)
        if pv is None:
            continue
        M[r], M[pv] = M[pv], M[r]
        inv = pow(M[r][col], -1, p)
        M[r] = [(x * inv) % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][col]:
                f = M[i][col]
                M[i] = [(x - f * y) % p for x, y in zip(M[i], M[r])]
        r += 1
    return r


def rango_Q(filas):
    from fractions import Fraction
    M = [[Fraction(x) for x in f] for f in filas]
    r = 0
    ncols = len(M[0]) if M else 0
    for col in range(ncols):
        pv = next((i for i in range(r, len(M)) if M[i][col] != 0), None)
        if pv is None:
            continue
        M[r], M[pv] = M[pv], M[r]
        for i in range(len(M)):
            if i != r and M[i][col] != 0:
                f = M[i][col] / M[r][col]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        r += 1
    return r


def matriz_momento(f, n):
    """[f(u^{-1} c)]_{u, c}: u en (Z/n)^x/+-1, c en Z/n."""
    filas = []
    for u in unidades_mod_pm(n):
        ui = pow(u, -1, n)
        filas.append([f[(ui * c) % n] for c in range(n)])
    return filas


def W1_predicho(exps, n, p):
    return rango_mod(matriz_momento(S_de(exps, n), n), p)


def sierra(n):
    return [0] + [2 * c - n for c in range(1, n)]


def sierra_desplazada(n):
    """hat s_n(c) = c (0<c<n/2), c-n (n/2<c<n), 0 en 0 y n/2."""
    r = [0] * n
    for c in range(1, n):
        if 2 * c < n:
            r[c] = c
        elif 2 * c > n:
            r[c] = c - n
    return r


def segmento_impar(m):
    return [2 * j - 1 for j in range(1, m + 1)]


def segmento_lleno(m):
    return list(range(1, m + 1))
