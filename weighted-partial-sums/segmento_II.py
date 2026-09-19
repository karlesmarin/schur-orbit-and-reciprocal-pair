#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""segmento_II.py -- vanishing weighted sums of completely multiplicative functions, computed.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

A single file, Python 3 standard library only.  No Sage, no numpy, no install.  Every number is an
exact integer (or an exact element of Z[zeta_N]); nothing is decided in floating point.

WHAT IT COMPUTES  (each routine is named after the statement of the note it implements)

    cero(m)             a completely multiplicative eps: N -> {+1,-1} with
                        sum_{n<=m} n eps(n) = 0, the sum recomputed exactly;
                        or, for m = 1, 2 mod 4, why none exists     (Thm. every admissible length)
    certificado(j)      re-verifies block j of the 104 interval certificates
                        that prove it for 5000 <= m < 10^8        (proof of the same theorem)
    digitos(q, s)       the digit set R_q of a Legendre symbol, the zeros of
                        S_chi below a bound, each checked by direct summation
                        AND by the digit recursion       (Thm. a digit recursion, Cor. zeros from digits)
    bloque(q, digs)     a loop of the recursion: the zeros q B (1+Q+...+Q^(t-1))       (Prop. blocks)
    semisuma(f)         the weighted half sum of EVERY odd character mod f
                        against its closed form, exactly in Z[zeta_N]     (Thm. the weighted half sum,
                                                      Cor. how many components vanish, Cor. Dirichlet)
    conrey(q)           the cells on which Conrey's sine series of the Legendre
                        symbol vanishes identically: they are R_q     (Sec. Conrey's sine series,
                                                      Prop. positivity empties the digit set)

USAGE
    python segmento_II.py autotest            recompute the worked examples of the note
    python segmento_II.py cero m [full]       a zero at length m (explicit up to CERO_MAX)
    python segmento_II.py certificado [j]     all 104 block certificates, or block j (1..104)
    python segmento_II.py digitos q [s] [X]   R_q and the zeros of S_chi up to X (s = +1 or -1)
    python segmento_II.py semisuma f          the half-sum formula for every odd character mod f
    python segmento_II.py conrey q            the null cells of F_q, q = 3 mod 4, q > 3

    As a library:  from segmento_II import cero, certificado, digitos, bloque, semisuma, conrey

WHAT IT DECLINES
    cero        above CERO_MAX it computes nothing and says so: existence there is the theorem
                (certificates to 10^8, the analytic tail beyond), not a computation of this tool.
    digitos     q not an odd prime; bounds above DIGITOS_MAX for direct summation.
    semisuma    f even, f < 3, or f above SEMISUMA_MAX.
    conrey      q not a prime = 3 mod 4 with q > 3: outside the identity the section uses.

Public domain (CC0).
"""
import itertools
import os
import sys
from fractions import Fraction
from math import gcd, isqrt

CERO_MAX = 10 ** 6       # explicit zeros are built and re-summed up to here (about a second)
DIGITOS_MAX = 10 ** 7    # direct summation bound for `digitos` (a few seconds)
SEMISUMA_MAX = 1000      # all odd characters mod f, exactly, up to here (seconds)
CERT_CSV = "II_certificados_A_bloques.csv"
CERT_TOP = 10 ** 8       # the certificates cover [5000, 10^8)


# --------------------------------------------------------------------------- basic arithmetic

def es_primo(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def factoriza(n):
    """{prime: exponent}."""
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


def orden_mult(a, n):
    a %= n
    k, x = 1, a
    while x != 1:
        x = x * a % n
        k += 1
    return k


def jacobi(a, n):
    """Jacobi symbol (a|n), n odd positive."""
    if n <= 0 or n % 2 == 0:
        raise ValueError("the Jacobi symbol needs n odd and positive")
    a %= n
    t = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                t = -t
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            t = -t
        a %= n
    return t if n == 1 else 0


def criba_spf(n):
    """smallest prime factor of every k <= n."""
    spf = list(range(n + 1))
    for i in range(2, isqrt(n) + 1):
        if spf[i] == i:
            for j in range(i * i, n + 1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


def primos_hasta(n):
    if n < 2:
        return []
    c = bytearray([1]) * (n + 1)
    c[0:2] = b"\x00\x00"
    for i in range(2, isqrt(n) + 1):
        if c[i]:
            c[i * i::i] = bytearray(len(range(i * i, n + 1, i)))
    return list(itertools.compress(range(n + 1), c))


def primos_intervalo(lo, hi, base=None):
    """primes in (lo, hi], by a segmented sieve."""
    if hi <= lo:
        return []
    if base is None:
        base = primos_hasta(isqrt(hi))
    ini = lo + 1
    c = bytearray([1]) * (hi - lo)
    for p in base:
        if p * p > hi:
            break
        s = max(p * p, (ini + p - 1) // p * p)
        if s <= hi:
            c[s - ini::p] = bytearray(len(range(s - ini, hi - lo, p)))
    if ini <= 1:
        for x in range(ini, 2):
            c[x - ini] = 0
    return list(itertools.compress(range(ini, hi + 1), c))


def numero_de_clases(D):
    """h(D), D < 0, counted as primitive reduced binary quadratic forms (a, b, c), b^2 - 4ac = D.
    Independent of every character sum below; used to cross-check them."""
    if D >= 0 or D % 4 not in (0, 1):
        raise ValueError("D must be a negative discriminant")
    h, a = 0, 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a):
                continue
            c = (b * b - D) // (4 * a)
            if c < a or (b < 0 and a == c):
                continue
            if gcd(gcd(a, abs(b)), c) != 1:
                continue
            h += 1
        a += 1
    return h


# ==================================================================== 1. every admissible length
#
#   S_eps(m) = sum_{n<=m} n eps(n).  Changing one sign changes it by an even amount, so
#   S_eps(m) = m(m+1)/2 (mod 2), and a zero needs m = 0, 3 (mod 4): the ADMISSIBLE lengths.

def admisible(m):
    return m % 4 in (0, 3)


def suma_ponderada(signos, m, spf=None):
    """S_eps(m) for the completely multiplicative eps with eps(p) = signos[p], exactly."""
    if spf is None:
        spf = criba_spf(max(m, 2))
    e = [0] * (m + 1)
    if m >= 1:
        e[1] = 1
    tot = 1 if m >= 1 else 0
    for c in range(2, m + 1):
        p = spf[c]
        e[c] = signos[p] * e[c // p]
        tot += c * e[c]
    return tot


def _mascaras(m, primos):
    """for each n <= m, the bitmask of the primes dividing n to an ODD power."""
    idx = {p: i for i, p in enumerate(primos)}
    om = [0] * (m + 1)
    for n in range(2, m + 1):
        x, msk = n, 0
        for p in primos:
            if p > x:
                break
            while x % p == 0:
                x //= p
                msk ^= 1 << idx[p]
        om[n] = msk
    return om


def soluciones_exhaustivas(m):
    """every sign pattern on the primes <= m with S_eps(m) = 0 (only for small m)."""
    if m > 60:
        raise ValueError("exhaustive search only for m <= 60")
    primos = primos_hasta(m)
    om = _mascaras(m, primos)
    out = []
    for msk in range(1 << len(primos)):
        s = 0
        for n in range(1, m + 1):
            s += -n if bin(msk & om[n]).count("1") % 2 else n
        if s == 0:
            out.append({p: (-1 if (msk >> i) & 1 else 1) for i, p in enumerate(primos)})
    return out


def _subconjunto_bitset(L, t):
    """a subset of L with sum t, by an exact bitset over all prefixes; None if there is none."""
    if t < 0 or t > sum(L):
        return None
    alc = [1]
    for p in L:
        alc.append(alc[-1] | (alc[-1] << p))
    if not (alc[-1] >> t) & 1:
        return None
    out = []
    for i in range(len(L) - 1, -1, -1):
        if not (alc[i] >> t) & 1:
            out.append(L[i])
            t -= L[i]
    return out if t == 0 else None


def _subconjunto(L, t, k=32):
    """a subset of the primes L of (m/2, m] with sum exactly t, or None.

    For long L: the largest primes are taken greedily while the remainder stays above the middle
    of the subset sums of k primes spread over L; those k then finish exactly with a bitset.  The
    answer is always re-checked by the caller, so a wrong subset cannot pass silently."""
    P = sum(L)
    if t < 0 or t > P:
        return None
    if len(L) <= 2 * k:
        return _subconjunto_bitset(L, t)
    paso = len(L) // k
    peq = L[::paso][:k]
    sp = set(peq)
    grandes = sorted((p for p in L if p not in sp), reverse=True)
    meta = sum(peq) // 2
    if t < meta:
        comp = _subconjunto(L, P - t, k)
        if comp is None:
            return None
        cs = set(comp)
        return [p for p in L if p not in cs]
    tomados, resto = [], t
    for p in grandes:
        if resto - p >= meta:
            tomados.append(p)
            resto -= p
    fin = _subconjunto_bitset(peq, resto)
    return None if fin is None else tomados + fin


# The starting signs at the primes <= sqrt(m), in the order they are tried: the note's route
# starts from eps(2) = eps(3) = -1 and changes at most two signs at primes <= 7.
VARIANTES = [{2: -1, 3: -1}, {2: -1, 3: -1, 5: -1}, {2: -1, 3: -1, 7: -1}, {2: -1}, {3: -1},
             {2: -1, 3: -1, 5: -1, 7: -1}, {2: -1, 5: -1}, {3: -1, 5: -1}]


def _ruta(m, spf, primos, base):
    """the constructive route of the proof: fixed signs at primes <= sqrt(m), the primes of
    (sqrt m, m/2] balanced greedily, and (RQ) solved EXACTLY with the primes of (m/2, m]."""
    raiz = isqrt(m)
    s = {p: base.get(p, 1) for p in primos if p <= raiz}
    e = [0] * (m + 1)
    e[1] = 1
    liso = [True] * (m + 1)      # all prime factors <= sqrt(m)
    for c in range(2, m + 1):
        p = spf[c]
        e[c] = (s[p] if p <= raiz else 1) * e[c // p]
        liso[c] = p <= raiz and liso[c // p]
    T = [0] * (raiz + 2)
    for k in range(1, raiz + 2):
        T[k] = T[k - 1] + (k * e[k] if k <= m else 0)
    R = sum(c * e[c] for c in range(1, m + 1) if liso[c])
    # A prime p of (sqrt m, m/2] meets n <= m only as n = p k with k <= m/p < sqrt(m)+1:
    # it moves the sum by eps(p) p T(m // p).
    medios = [p for p in primos if raiz < p <= m // 2]
    for p in sorted(medios, key=lambda p: -abs(p * T[m // p])):
        w = p * T[m // p]
        if abs(R - w) < abs(R + w):
            R -= w
            s[p] = -1
        else:
            R += w
            s[p] = 1
    # (RQ): S = R + sum_{q in (m/2, m]} q eps(q).  The primes with eps = -1 must sum to (P + R)/2.
    L = [p for p in primos if m // 2 < p <= m]
    P = sum(L)
    if (P + R) % 2 or abs(R) > P:
        return None
    neg = _subconjunto(L, (P + R) // 2)
    if neg is None:
        return None
    ns = set(neg)
    for p in L:
        s[p] = -1 if p in ns else 1
    return s


def cero(m, limite=None):
    """Theorem 'every admissible length', made explicit at one length m.

    Returns (signos, S, how): signos = {p: eps(p)} on every prime p <= m, S = the weighted sum
    RECOMPUTED from those signs (always 0 when returned), and how it was found.  For m = 1, 2
    mod 4 returns (None, T_m, reason): no eps exists, by parity.  Above `limite` (CERO_MAX by
    default) it raises ValueError: it does not claim what it has not computed."""
    limite = CERO_MAX if limite is None else limite
    if m < 1:
        raise ValueError("m must be a positive integer")
    if not admisible(m):
        T = m * (m + 1) // 2
        return None, T, ("no eps exists: every eps gives S_eps(m) = m(m+1)/2 = %d (mod 2), which "
                         "is odd, since m = %d (mod 4)" % (T, m % 4))
    if m > limite:
        # GUARD: above-range
        raise ValueError(
            "m = %d is above the explicit range (%d): nothing is computed here.  A zero EXISTS by "
            "the theorem 'every admissible length' -- interval certificates for 5000 <= m < 10^8 "
            "(`certificado`) and the analytic tail for m >= 10^8 -- but this tool does not "
            "exhibit it." % (m, limite))
    spf = criba_spf(m)
    primos = [p for p in range(2, m + 1) if spf[p] == p]
    how = None
    if m <= 40:
        sols = soluciones_exhaustivas(m)
        if sols:
            s, how = sols[0], "exhaustive search over the %d sign patterns" % (2 ** len(primos))
    else:
        for i, b in enumerate(VARIANTES):
            s = _ruta(m, spf, primos, b)
            if s is not None:
                how = "route of the proof, starting signs %s" % (
                    ", ".join("eps(%d)=-1" % p for p in sorted(b)))
                break
    if how is None:
        raise RuntimeError("m = %d: the constructive route found no witness (existence is the "
                           "theorem; this tool did not find one)" % m)
    S = suma_ponderada(s, m, spf)
    if S != 0:
        raise AssertionError("internal error: the witness for m = %d sums to %d" % (m, S))
    return s, S, how


# ------------------------------------------------------------------- the 104 block certificates

def bloques_A():
    """[(a_j, b_j)]: a_1 = 5000, b_j = min(floor(11 a_j / 10), 10^8 - 1), a_{j+1} = b_j + 1."""
    out, a = [], 5000
    while a <= CERT_TOP - 1:
        b = min(11 * a // 10, CERT_TOP - 1)
        out.append((a, b))
        a = b + 1
    return out


def _lee_csv(ruta=None):
    ruta = ruta or os.path.join(os.path.dirname(os.path.abspath(__file__)), CERT_CSV)
    if not os.path.exists(ruta):
        raise ValueError("%s not found next to segmento_II.py: the certificates travel with it"
                         % CERT_CSV)
    filas = []
    with open(ruta, encoding="utf-8") as fh:
        cab = fh.readline().strip().split(",")
        if cab != ["a", "b", "N", "pares", "L", "D", "cota_sobre_D2"]:
            raise ValueError("unexpected header in %s: %s" % (CERT_CSV, cab))
        for linea in fh:
            if linea.strip():
                x = linea.strip().split(",")
                filas.append(tuple(int(v) for v in x[:6]) + (x[6],))
    return filas


def _log2_techo(b):
    """ceil(log2 b), b >= 2, in integers.  Then log b < 7e/10, since 0.7 > log 2."""
    return (b - 1).bit_length()


def certificado(j, filas=None, base=None):
    """Re-verify block j (1..104) of the proof of 'every admissible length' in [5000, 10^8).

    Q = primes of (b/2, a]; a prime q0 and disjoint pairs (u, v) of Q with d = (v-u)/2,
    1 <= d_i <= 1 + d_1 + ... + d_{i-1}, H = sum d, U = sum u, 2H - q0 + 2 >= b (Lemma 'intervals
    of subset sums' (a), (b)); alpha = U + q0 - 1, D = P_Q - 2 alpha; and the integer inequality
    90 D^2 > b^3 (130 + 7e), e = ceil(log2 b), which is B_*(b) < D.  The pairs are rebuilt by the
    same greedy rule and every number is compared with the archived row."""
    filas = _lee_csv() if filas is None else filas
    bl = bloques_A()
    if len(filas) != len(bl):
        raise ValueError("the CSV has %d rows and the block chain has %d" % (len(filas), len(bl)))
    if not isinstance(j, int) or not 1 <= j <= len(bl):
        # GUARD: block-index
        raise ValueError("block index must be 1..%d, got %r" % (len(bl), j))
    a, b = bl[j - 1]
    ca, cb, cN, cpares, cL, cD, ccota = filas[j - 1]
    fallos = []
    if (ca, cb) != (a, b):
        fallos.append("block ends (%d, %d) in the CSV, (%d, %d) by the rule" % (ca, cb, a, b))
    Q = primos_intervalo(b // 2, a, base)
    Qs = set(Q)
    PQ = sum(Q)
    q0 = Q[0]
    usados = {q0}
    H = U = 0
    pares = []
    while 2 * H - q0 + 2 < b:
        hecho = False
        for d in range(min(H + 1, (a - Q[0]) // 2), 0, -1):
            for u in Q:
                v = u + 2 * d
                if v > a:
                    break
                if u in usados or v in usados or v not in Qs:
                    continue
                pares.append((u, v))
                usados.update((u, v))
                H += d
                U += u
                hecho = True
                break
            if hecho:
                break
        if not hecho:
            fallos.append("no admissible pair left at H = %d" % H)
            break
    # the lemma's hypotheses, one by one
    Hc = 0
    for (u, v) in pares:
        dj = (v - u) // 2
        if not ((v - u) % 2 == 0 and 1 <= dj <= Hc + 1 and b // 2 < u < v <= a
                and u in Qs and v in Qs):
            fallos.append("pair (%d, %d) breaks the pairs lemma" % (u, v))
        Hc += dj
    if not (Hc == H and 2 * H - q0 + 2 >= b):
        fallos.append("2H - q0 + 2 >= b fails")
    alpha = U + q0 - 1
    D = PQ - 2 * alpha
    e = _log2_techo(b)
    ok = 90 * D * D > b ** 3 * (130 + 7 * e)
    if not ok:
        fallos.append("90 D^2 > b^3 (130 + 7e) FAILS")
    # control: the certificate is not vacuous -- without its last pair the condition fails
    if pares:
        u, v = pares[-1]
        if 2 * (H - (v - u) // 2) - q0 + 2 >= b:
            fallos.append("the last pair is not needed: the control does not bite")
    cota = Fraction(b ** 3 * (130 + 7 * e), 90 * D * D)
    for nombre, arch, calc in (("N", cN, len(Q)), ("pairs", cpares, len(pares)),
                               ("L", cL, alpha), ("D", cD, D)):
        if arch != calc:
            fallos.append("%s = %d in the CSV, %d recomputed" % (nombre, arch, calc))
    if ccota != "%.6e" % float(cota):
        fallos.append("bound/D^2 = %s in the CSV, %.6e recomputed" % (ccota, float(cota)))
    return {"j": j, "a": a, "b": b, "N": len(Q), "q0": q0, "pairs": len(pares), "H": H,
            "alpha": alpha, "D": D, "e": e, "ratio": cota,
            "prefix16": 90 * D * D > 16 * b ** 3 * (130 + 7 * e), "failures": fallos}


def certificados(indices=None, eco=None):
    """re-verify several blocks (all by default); returns (results, total failures)."""
    filas = _lee_csv()
    base = primos_hasta(isqrt(CERT_TOP) + 1)
    res = []
    for j in (indices or range(1, len(filas) + 1)):
        r = certificado(j, filas, base)
        res.append(r)
        if eco:
            eco(r)
    return res, sum(len(r["failures"]) for r in res)


# =========================================================== 2. along one function: characters
#
#   chi = chi_{q,s}: the Legendre symbol (.|q) on integers prime to q, chi(q) = s = +-1, extended
#   completely multiplicatively.  C(r) = sum_{r'<=r} chi(r'), T(r) = sum_{r'<=r} r' chi(r'),
#   T_q = T(q-1), for 0 <= r < q.

def tablas(q):
    """(L, C, T): Legendre table mod q and its unweighted and weighted partial sums."""
    L = [0] + [-1] * (q - 1)
    for r in range(1, q):
        L[r * r % q] = 1
    C, T, c, t = [0] * q, [0] * q, 0, 0
    for r in range(q):
        c += L[r]
        t += r * L[r]
        C[r], T[r] = c, t
    return L, C, T


def _exige_primo_impar(q):
    if not (isinstance(q, int) and q >= 3 and es_primo(q)):
        # GUARD: odd-prime
        raise ValueError("q = %r must be an odd prime (the Legendre symbol mod q)" % (q,))


def conjunto_R(q):
    """R_q = {0 <= r < q : T_q + q C(r) = 0 = T(r)}  (Cor. zeros from digits)."""
    _exige_primo_impar(q)
    L, C, T = tablas(q)
    return [r for r in range(q) if T[q - 1] + q * C[r] == 0 and T[r] == 0]


def S_directa(q, s, m):
    """sum_{n<=m} n chi_{q,s}(n), term by term."""
    L = tablas(q)[0]
    tot = 0
    for n in range(1, m + 1):
        x, sg = n, 1
        while x % q == 0:
            x //= q
            sg *= s
        tot += n * sg * L[x % q]
    return tot


def S_recursion(q, s, m, tab=None):
    """the same sum by Thm. 'a digit recursion', read on the base-q digits of m:
    S(qk + r) = s q S(k) + k (T_q + q C(r)) + T(r)."""
    L, C, T = tab or tablas(q)
    dig = []
    x = m
    while x:
        dig.append(x % q)
        x //= q
    S = k = 0
    for r in reversed(dig):
        S = s * q * S + k * (T[q - 1] + q * C[r]) + T[r]
        k = q * k + r
    return S


def clase_q(q):
    """h(-q) for a prime q = 3 mod 4, q > 3, by three routes that must agree: reduced forms, the
    half sum of Jacobi and Dirichlet h = sum_{n<q/2} (n|q) / (2 - (2|q)), and T_q = -q h."""
    L, C, T = tablas(q)
    h = numero_de_clases(-q)
    num = C[(q - 1) // 2]
    den = 2 - L[2]
    jd = Fraction(num, den)
    if jd != h or T[q - 1] != -q * h:
        raise AssertionError("class number routes disagree at q = %d: forms %d, Jacobi-Dirichlet "
                             "%s, -T_q/q = %s" % (q, h, jd, Fraction(-T[q - 1], q)))
    return h


def arbol_digitos(q, R, X):
    """every m in [1, X] whose base-q digits all lie in R."""
    out, frente = [], [0]
    while frente:
        nuevo = []
        for k in frente:
            for r in R:
                m = q * k + r
                if 0 < m <= X:
                    out.append(m)
                    nuevo.append(m)
        frente = nuevo
    return sorted(set(out))


def ceros_hasta(q, s, X):
    """(zeros of S_chi in [1, X] by direct summation, number of m <= X at which the digit
    recursion was checked against the direct sum, failures of that check)."""
    if X > DIGITOS_MAX:
        # GUARD: digitos-bound
        raise ValueError("X = %d is above the direct-summation bound %d; use S_recursion for "
                         "single large m" % (X, DIGITOS_MAX))
    L, C, T = tablas(q)
    Tq = T[q - 1]
    guarda = X // q + 1
    Sk = [0] * (guarda + 1)
    tot, ceros, malos = 0, [], 0
    for n in range(1, X + 1):
        x, sg = n, 1
        while x % q == 0:
            x //= q
            sg *= s
        tot += n * sg * L[x % q]
        if n <= guarda:
            Sk[n] = tot
        k, r = divmod(n, q)
        if s * q * Sk[k] + k * (Tq + q * C[r]) + T[r] != tot:
            malos += 1
        if tot == 0:
            ceros.append(n)
    return ceros, X, malos


def digitos(q, s=1, X=10 ** 6):
    """Thm. 'a digit recursion' and Cor. 'zeros from digits' at one prime q and sign s.

    Returns a dict with R_q, h (q = 3 mod 4, q > 3), the zeros below X (each by direct summation
    and by the recursion), and every part of the corollary that applies, checked."""
    _exige_primo_impar(q)
    if s not in (1, -1):
        raise ValueError("s = chi(q) must be +1 or -1")
    tab = tablas(q)
    L, C, T = tab
    R = conjunto_R(q)
    out = {"q": q, "s": s, "X": X, "R": R, "T_q": T[q - 1], "partes": [], "fallos": []}
    h = clase_q(q) if (q % 4 == 3 and q > 3) else None
    out["h"] = h
    ceros, nchk, malos = ceros_hasta(q, s, X)
    out["ceros"] = ceros
    if malos:
        out["fallos"].append("the digit recursion differs from the direct sum at %d m <= X" % malos)
    out["recursion_checked_at"] = nchk
    zs = set(ceros)
    no_rec = [z for z in ceros if S_recursion(q, s, z, tab) != 0]
    if no_rec:
        out["fallos"].append("zeros by direct summation that the recursion does not confirm: %s"
                             % no_rec[:5])
    # the digit tree: every m with digits in R_q is a zero, for both s
    arbol = arbol_digitos(q, R, X)
    fuera = [m for m in arbol if m not in zs]
    if fuera:
        out["fallos"].append("digit-tree numbers that are not zeros: %s" % fuera[:5])
    out["arbol"] = arbol
    if q % 4 == 1:
        ok = 0 in R and q - 1 in R
        out["partes"].append(("(a) q = 1 mod 4: 0 and q-1 lie in R_q", ok))
        if not ok:
            out["fallos"].append("(a) fails")
    if q % 8 == 7:
        ok = (q - 1) // 2 in R
        zj = [(q ** j - 1) // 2 for j in range(1, 13)]
        rec = all(S_recursion(q, s, z, tab) == 0 for z in zj)
        dirz = all(z in zs for z in zj if z <= X)
        out["partes"].append(("(b) q = 7 mod 8: (q-1)/2 in R_q; S((q^j-1)/2) = 0 for j <= 12 by the "
                              "recursion, and by direct summation for the %d of them <= X"
                              % sum(1 for z in zj if z <= X), ok and rec and dirz))
        if not (ok and rec and dirz):
            out["fallos"].append("(b) fails")
    if s == 1 and ((q % 4 == 3 and q > 3 and h == 1) or q == 3):
        base = q if q > 3 else 27
        cad = [base]
        while len(cad) < 12:
            cad.append(base * (cad[-1] + 1))
        rec = all(S_recursion(q, s, z, tab) == 0 for z in cad)
        dirz = all(z in zs for z in cad if z <= X)
        out["partes"].append(("(c) %s: the chain z_0 = %d, z_(i+1) = %d(z_i + 1) -- 12 links by the "
                              "recursion, %d by direct summation"
                              % ("h(-q) = 1, s = 1" if q > 3 else "q = 3, s = 1", base, base,
                                 sum(1 for z in cad if z <= X)), rec and dirz))
        if not (rec and dirz):
            out["fallos"].append("(c) fails")
    if s == -1:
        pares = [z for z in ceros if z * q * q <= X]
        ok = all(z * q * q in zs for z in pares)
        out["partes"].append(("(e) s = -1: q^2 z is a zero for each of the %d zeros z <= X/q^2"
                              % len(pares), ok))
        if not ok:
            out["fallos"].append("(e) fails")
    if q % 8 == 3 and q > 3:
        out["partes"].append(("q = 3 mod 8: R_q is %s (the note measures it empty for every "
                              "q < 10^5 and does not prove it)" % ("empty" if not R else
                                                                   "NOT empty"), not R))
    return out


def bloque(q, digs, t_max=3, h=None):
    """Prop. 'blocks': q = 3 mod 4, q > 3, s = 1, D(r) = C(r) - h.  If B has base-q digits digs
    (most significant first) with sum D(d_i) = 0 and S(B) = h B, every m_t = q B (1 + Q + ... +
    Q^(t-1)), Q = q^len(digs), is a zero.  Returns (hypotheses hold, [(m_t, S(m_t))])."""
    _exige_primo_impar(q)
    if q % 4 != 3 or q == 3:
        # GUARD: bloque-hyp
        raise ValueError("the block proposition needs q = 3 mod 4, q > 3")
    tab = tablas(q)
    L, C, T = tab
    h = clase_q(q) if h is None else h
    B = 0
    for d in digs:
        B = B * q + d
    hip = sum(C[d] - h for d in digs) == 0 and S_recursion(q, 1, B, tab) == h * B
    Qb = q ** len(digs)
    zs = []
    for t in range(1, t_max + 1):
        m = q * B * sum(Qb ** i for i in range(t))
        zs.append((m, S_recursion(q, 1, m, tab)))
    return hip, B, zs


# ============================================================ 3. the weighted half sum, exactly
#
#   For odd f and an odd Dirichlet character chi mod f,
#       sum_{c=1}^{(f-1)/2} c chi(c) = (f/2) B_{1,chi} (chi(2)^{-1} - 1),  f B_{1,chi} = sum a chi(a).
#   Multiplied by 2 chi(2) it is an identity in Z[zeta_N] with no denominator:
#       2 chi(2) T = A (1 - chi(2)),  A = sum_{a=1}^{f-1} a chi(a).

def ciclotomico(n):
    """Phi_n over Z, low degree first."""
    num = [-1] + [0] * (n - 1) + [1]
    for d in range(1, n):
        if n % d == 0:
            num = _div_exacta(num, ciclotomico(d))
    return num


def _div_exacta(a, b):
    a = a[:]
    q = [0] * (len(a) - len(b) + 1)
    for d in range(len(a) - len(b), -1, -1):
        c = a[d + len(b) - 1]
        q[d] = c
        if c:
            for i, y in enumerate(b):
                a[i + d] -= c * y
    return q


_PHI = {}


def es_cero_ciclotomico(coef, N):
    """is sum_k coef[k] zeta_N^k = 0 ?  (coef has length N; exact reduction modulo Phi_N)."""
    if N not in _PHI:
        _PHI[N] = ciclotomico(N)
    phi = _PHI[N]
    g = len(phi) - 1
    a = list(coef)
    for d in range(len(a) - 1, g - 1, -1):
        c = a[d]
        if c:
            base = d - g
            for i, y in enumerate(phi):
                if y:
                    a[base + i] -= c * y
    return not any(a[:g])


class Caracteres(object):
    """the Dirichlet characters mod an odd f, as exponent vectors on a product of cyclic groups.
    chi_a(x) = zeta_N^{e_a(x)}, N the exponent of (Z/f)^x, e_a(x) = sum a_i log_i(x) N/n_i."""

    def __init__(self, f):
        self.f = f
        self.comp = []
        for p, k in sorted(factoriza(f).items()):
            m = p ** k
            n = (p - 1) * p ** (k - 1)
            g = next(g for g in range(2, m + 1) if g % p and
                     all(pow(g, n // r, m) != 1 for r in factoriza(n)))
            logt, x = {}, 1
            for i in range(n):
                logt[x] = i
                x = x * g % m
            self.comp.append((p, m, n, logt))
        self.N = 1
        for _, _, n, _ in self.comp:
            self.N = self.N * n // gcd(self.N, n)
        self.unidades = [x for x in range(1, f) if gcd(x, f) == 1]
        self.logs = {x: [lt[x % m] for (_, m, _, lt) in self.comp] for x in self.unidades}
        # for primitivity: the units = 1 mod f/p, for each prime p | f
        self.nucleos = []
        for p, _, _, _ in self.comp:
            g = f // p
            self.nucleos.append([x for x in self.unidades if x % g == 1 % g])

    def todos(self):
        return itertools.product(*[range(n) for _, _, n, _ in self.comp])

    def e(self, a, x):
        """the exponent of chi_a(x) in Z/N, or None when gcd(x, f) > 1."""
        lg = self.logs.get(x % self.f)
        if lg is None:
            return None
        return sum(ai * li * (self.N // n) for ai, li, (_, _, n, _) in zip(a, lg, self.comp)) % self.N

    def impar(self, a):
        return self.e(a, self.f - 1) == self.N // 2

    def primitivo(self, a):
        return all(any(self.e(a, x) != 0 for x in nuc) for nuc in self.nucleos)

    def trivial(self, a):
        return all(ai == 0 for ai in a)

    def suma(self, a, m):
        """coefficients of sum_{c=1}^{m} c chi_a(c) on the powers of zeta_N."""
        v = [0] * self.N
        for c in range(1, m + 1):
            k = self.e(a, c)
            if k is not None:
                v[k] += c
        return v


def semisuma(f, limite=None, fuera_diagonal=False, pares=False):
    """Thm. 'the weighted half sum' for EVERY odd character mod f, exactly in Z[zeta_N].

    Returns counts: odd characters, closed form verified, primitive ones, criterion 'T = 0 iff
    chi(2) = 1' right on the primitive ones and overall, and the number of vanishing T.  With
    fuera_diagonal=True the criterion is tested one step off the diagonal, m = (f-1)/2 - 1, and
    with pares=True on the even non-trivial characters: those are DECOYS, which must fail."""
    limite = SEMISUMA_MAX if limite is None else limite
    if not isinstance(f, int) or f < 3 or f % 2 == 0:
        # GUARD: semisuma-odd
        raise ValueError("f = %r: the theorem is for odd f >= 3" % (f,))
    if f > limite:
        # GUARD: semisuma-range
        raise ValueError("f = %d is above %d: the exact sweep over all characters is not run here"
                         % (f, limite))
    G = Caracteres(f)
    N = G.N
    m = (f - 1) // 2 - (1 if fuera_diagonal else 0)
    r = {"f": f, "N": N, "caracteres": 0, "forma_ok": 0, "primitivos": 0, "crit_prim_ok": 0,
         "crit_ok": 0, "ceros": 0, "ceros_chi2_1": 0}
    for a in G.todos():
        if pares:
            if G.impar(a) or G.trivial(a):
                continue
        elif not G.impar(a):
            continue
        if m < 1:
            continue
        r["caracteres"] += 1
        T = G.suma(a, m)
        k2 = G.e(a, 2)
        cero_T = es_cero_ciclotomico(T, N)
        r["ceros"] += cero_T
        r["ceros_chi2_1"] += cero_T and k2 == 0
        crit = cero_T == (k2 == 0)
        r["crit_ok"] += crit
        if not fuera_diagonal and not pares:
            A = G.suma(a, f - 1)
            # 2 chi(2) T - A (1 - chi(2)) = 0 ?
            d = [0] * N
            for k in range(N):
                d[(k + k2) % N] += 2 * T[k] + A[k]
                d[k] -= A[k]
            r["forma_ok"] += es_cero_ciclotomico(d, N)
            if G.primitivo(a):
                r["primitivos"] += 1
                r["crit_prim_ok"] += crit
    return r


def cuantos(f):
    """Cor. 'how many components vanish': for f PRIME, 0 if -1 is in <2>, else d/t.  Declines on
    composite f, where imprimitive characters add zeros the count does not see (first at f = 21)."""
    if not (isinstance(f, int) and f >= 5 and es_primo(f)):
        # GUARD: cuantos-prime
        raise ValueError("f = %r: the count d/t is proved for PRIME f >= 5 only (at f = 21 it "
                         "predicts 1 and there are 2)" % (f,))
    d = (f - 1) // 2
    t = orden_mult(2, f)
    en = (f - 1) in {pow(2, j, f) for j in range(t)}
    return 0 if en else d // t


def dirichlet(f):
    """Cor. 'Dirichlet's family': f = 3 mod 4 squarefree, chi = (.|f) (Jacobi), odd and of
    conductor f.  Returns (sum_{c<=(f-1)/2} c chi(c), the predicted value): f h(-f) 2/w(-f) when
    f = 3 mod 8, and 0 when f = 7 mod 8; h(-f) by reduced forms."""
    if f % 4 != 3 or any(k > 1 for k in factoriza(f).values()):
        # GUARD: dirichlet-hyp
        raise ValueError("f = %d: the corollary needs f = 3 mod 4 squarefree" % f)
    T = sum(c * jacobi(c, f) for c in range(1, (f - 1) // 2 + 1))
    if f % 8 == 7:
        return T, 0
    w = 6 if f == 3 else 2
    return T, Fraction(f * numero_de_clases(-f) * 2, w)


# ================================================================== 4. Conrey's sine series
#
#   For a prime q = 3 mod 4, q > 3, h = h(-q):  f_q(x) = (2 pi^2 / q^{3/2}) F_q(qx) with
#   F_q(y) = T(r) + y (h - C(r)), r = floor(y).  F_q is linear on each cell [r, r+1] and vanishes
#   identically on it exactly when C(r) = h and T(r) = 0.  Only F_q is computed here, exactly; the
#   series itself is not summed.

def conrey(q):
    """the null cells of F_q, compared with R_q, and whether Prop. 'positivity empties the digit
    set' (q = 3 mod 8, q > 3) applies at q."""
    _exige_primo_impar(q)
    if q % 4 != 3 or q == 3:
        # GUARD: conrey-hyp
        raise ValueError("q = %d: the identity with F_q is used for primes q = 3 mod 4, q > 3"
                         % q)
    L, C, T = tablas(q)
    h = clase_q(q)
    nulas = []
    for r in range(q):
        izq = T[r] + r * (h - C[r])          # F_q at y = r, from the cell of r
        der = T[r] + (r + 1) * (h - C[r])    # F_q at y = r + 1, from the cell of r
        if izq == 0 and der == 0:
            nulas.append(r)
    # continuity of F_q across every node: the formula of cell r and of cell r+1 agree at r+1
    continua = all(T[r] + (r + 1) * (h - C[r]) == T[r + 1] + (r + 1) * (h - C[r + 1])
                   for r in range(q - 1))
    # the mirror used in the proof of the proposition
    espejo = all(C[q - 1 - r] == C[r] and T[q - 1 - r] == q * (C[r] - h) - T[r]
                 for r in range(q))
    R = conjunto_R(q)
    return {"q": q, "h": h, "nulas": nulas, "R": R, "igual": nulas == R, "continua": continua,
            "espejo": espejo, "prop_aplica": q % 8 == 3,
            "celdas_x": [(Fraction(r, q), Fraction(r + 1, q)) for r in nulas]}


# ------------------------------------------------------------------------------------ autotest

def _linea(ok, texto):
    print("   %-78s %s" % (texto, "ok" if ok else "*** FAILS ***"))
    return 0 if ok else 1


def _declina(fn, *args):
    try:
        fn(*args)
    except ValueError:
        return True
    return False


def autotest():
    f = 0

    print("every admissible length: explicit zeros, the sum recomputed")
    malos = [m for m in range(3, 401) if admisible(m) and cero(m)[1] != 0]
    f += _linea(not malos, "all %d admissible m <= 400 carry an explicit zero"
                % sum(1 for m in range(3, 401) if admisible(m)))
    for m in (12, 15, 19, 23):
        f += _linea(cero(m)[2].startswith("exhaustive"), "m = %d: by exhaustive search" % m)
    for m in (5000, 20000):
        s, S, how = cero(m)
        f += _linea(S == 0 and how.startswith("route"), "m = %d: S = %d, %s" % (m, S, how))
    ams = [m for m in range(3, 25) if admisible(m)]
    a = [len(soluciones_exhaustivas(m)) for m in ams]
    f += _linea(a == [1, 1, 1, 1, 2, 2, 2, 4, 1, 9, 5, 14],
                "a(m) on admissible m <= 24: %s" % a)
    print("parity decoys: they must find nothing, or fail")
    inad = [m for m in range(1, 23) if not admisible(m)]
    f += _linea(all(not soluciones_exhaustivas(m) for m in inad),
                "no zero at the %d inadmissible m <= 22 (exhaustive search)" % len(inad))
    f += _linea(all(cero(m)[0] is None for m in inad), "cero() answers 'none, by parity' at them")
    ocup = {m: admisible(m) and cero(m)[1] == 0 for m in range(2, 41)}
    j2 = sum(ocup[m] == (jacobi(2, 2 * m + 1) == 1) for m in range(2, 41))
    j3 = sum(ocup[m] == (jacobi(3, 2 * m + 1) == 1) for m in range(2, 41))
    f += _linea(j2 == 39 and j3 == 19,
                "rows 2..40 occupied iff (2|n)=1: %d of 39; decoy (3|n): %d of 39" % (j2, j3))
    vals35 = [jacobi(c, 35) for c in range(1, 7)]
    t35 = sum(c * v for c, v in enumerate(vals35, 1))
    f += _linea(t35 == 0 and jacobi(2, 13) == -1 and 0 in vals35,
                "f=35, m=6: T = 0 with (2|13) = -1, because chi(5) = 0 (not a +-1 function)")
    f += _linea(_declina(cero, CERO_MAX + 3), "above CERO_MAX: declines, computes nothing")
    f += _linea(_declina(cero, 0), "m = 0: declines")

    print("the block certificates of [5000, 10^8)")
    res, nf = certificados([1, 2, 3, 104])
    for r in res:
        f += _linea(not r["failures"], "block %3d [%d, %d]: %d pairs, B*^2/D^2 <= %.4f" % (
            r["j"], r["a"], r["b"], r["pairs"], float(r["ratio"])))
    f += _linea(float(res[0]["ratio"]) <= 0.445, "the first block's bound is at most 0.445")
    f += _linea(len(bloques_A()) == 104 and bloques_A()[-1][1] == CERT_TOP - 1,
                "the rule gives 104 blocks ending at 10^8 - 1")
    f += _linea(_declina(certificado, 0) and _declina(certificado, 105),
                "block indices 0 and 105: decline")

    print("the digit recursion and the digit sets R_q")
    esper = {5: [0, 4], 13: [0, 12], 29: [0, 4, 24, 28], 7: [3], 23: [11], 103: [47, 51, 55],
             11: [], 19: [], 43: [], 67: [], 163: [], 59: []}
    for q, R in sorted(esper.items()):
        f += _linea(conjunto_R(q) == R, "R_%d = %s" % (q, conjunto_R(q)))
    for q, ej in ((5, [4, 20, 24, 100]), (13, [12, 156, 168]), (29, [4, 24, 28, 116]),
                  (7, [3, 24, 171, 7, 56, 399]), (11, [11, 132, 1463, 16104]),
                  (131, [27, 5660, 743583])):
        tab = tablas(q)
        ok = all(S_directa(q, 1, m) == 0 and S_recursion(q, 1, m, tab) == 0 for m in ej)
        f += _linea(ok, "q=%d, s=1: %s are zeros (direct and recursion)" % (q, ej))
    for q in (5, 13, 29):
        d1, d2 = digitos(q, 1, 10 ** 6), digitos(q, -1, 10 ** 6)
        ok = (d1["ceros"] == d1["arbol"] == d2["ceros"] and not d1["fallos"] and not d2["fallos"])
        f += _linea(ok, "q=%d: the zeros <= 10^6 are the digit tree, %d of them, both s"
                    % (q, len(d1["ceros"])))
    d = digitos(23, 1, 10 ** 6)
    f += _linea(d["ceros"] == [11, 264, 6083, 139920] and not d["fallos"],
                "q=23, s=1: zeros <= 10^6 are %s = (23^j-1)/2" % d["ceros"])
    d = digitos(59, 1, 4600000)
    f += _linea(d["ceros"] == [4591439] and not d["fallos"],
                "q=59, s=1: the only zero <= 4.6*10^6 is %s" % d["ceros"])
    for q in (7, 23, 31, 47, 71, 79, 103):
        tab = tablas(q)
        zj = [(q ** j - 1) // 2 for j in range(1, 9)]
        ok = all(S_recursion(q, s, z, tab) == 0 for z in zj for s in (1, -1))
        ok = ok and all(S_directa(q, s, z) == 0 for z in zj if z <= 3 * 10 ** 5 for s in (1, -1))
        f += _linea(ok, "q=%d = 7 mod 8: S((q^j-1)/2) = 0, j <= 8, both s" % q)
    ok = all(S_recursion(q, -1, q * q * m) == q * q * S_directa(q, -1, m)
             for q in (11, 19, 23) for m in range(1, 3001, 37))
    f += _linea(ok, "(e) s=-1: S(q^2 m) = q^2 S(m) in the 246 cases q=11,19,23, m=1..3000 step 37")
    for q in (7, 11, 19, 43, 67, 163, 3):
        dq = digitos(q, 1, 3 * 10 ** 5)
        f += _linea(not dq["fallos"] and all(ok for _, ok in dq["partes"]),
                    "q=%d, s=1: every part of the corollary that applies holds, recursion ok "
                    "at every m <= 3*10^5" % q)
    q, tab = 131, tablas(131)
    fam = [27 + 43 * sum(q ** j for j in range(1, t + 1)) for t in range(6)]
    f += _linea(all(S_recursion(q, 1, m, tab) == 0 for m in fam),
                "q=131: 27 + 43(131 + ... + 131^t) is a zero for t <= 5 (recursion)")
    hip, B, zs = bloque(59, (1, 22, 21, 0))
    f += _linea(hip and B == 283200 and zs[0][0] == 16708800 and all(z == 0 for _, z in zs),
                "block q=59, B=(1,22,21,0)=283200: S(B)=3B, zeros %s" % [m for m, _ in zs])
    f += _linea(S_directa(59, 1, 4591439) == 0 and S_recursion(59, 1, 283200) == 849600,
                "q=59: S(4591439) = 0 and S(283200) = 849600 by direct summation / recursion")
    hip, B, zs = bloque(347, (40, 49, 226, 119, 164, 210, 117, 28, 302, 58, 264, 0, 0), 2)
    f += _linea(hip and all(z == 0 for _, z in zs), "block q=347, thirteen digits, h=5: zeros")
    bad = [q for q in primos_hasta(2000) if q > 3 and q % 4 == 3 and
           (clase_q(q) is None or (q % 8 == 3 and conjunto_R(q)) or
            (q % 8 == 7 and (q - 1) // 2 not in conjunto_R(q)) or
            (q % 8 == 3 and sum(c * tablas(q)[0][c] for c in range(1, (q + 1) // 2))
             != q * clase_q(q)))]
    f += _linea(not bad, "q = 3 mod 4, 3 < q < 2000: h by 3 routes; R_q empty at 3 mod 8, "
                         "(q-1)/2 in R_q at 7 mod 8; half sum = q h at 3 mod 8")
    f += _linea(_declina(digitos, 9) and _declina(digitos, 2) and _declina(digitos, 7, 1,
                                                                          DIGITOS_MAX + 1),
                "q = 9, q = 2, X above the bound: decline")
    f += _linea(_declina(bloque, 13, (1,)), "block at q = 13 (1 mod 4): declines")

    print("the weighted half sum, every odd character, exactly in Z[zeta_N]")
    tot = {"caracteres": 0, "forma_ok": 0, "primitivos": 0, "crit_prim_ok": 0, "crit_ok": 0}
    for ff in range(3, 122, 2):
        r = semisuma(ff)
        for k in tot:
            tot[k] += r[k]
    f += _linea(tot["caracteres"] == 1510 and tot["forma_ok"] == 1510,
                "odd f <= 121: %d odd characters, closed form in %d" % (tot["caracteres"],
                                                                      tot["forma_ok"]))
    f += _linea(tot["primitivos"] == 1236 and tot["crit_prim_ok"] == 1236,
                "criterion T = 0 iff chi(2) = 1 on the %d primitive ones: %d right"
                % (tot["primitivos"], tot["crit_prim_ok"]))
    f += _linea(tot["caracteres"] - tot["crit_ok"] == 19,
                "the same criterion without primitivity fails %d times (imprimitive, B_1 = 0)"
                % (tot["caracteres"] - tot["crit_ok"]))
    s1 = [semisuma(ff, fuera_diagonal=True) for ff in range(5, 60, 2)]
    s2 = [semisuma(ff, pares=True) for ff in range(5, 60, 2)]
    n1 = sum(r["caracteres"] - r["crit_ok"] for r in s1)
    n2 = sum(r["caracteres"] - r["crit_ok"] for r in s2)
    f += _linea(n1 == 15 and n2 == 9, "decoys, f < 60: one step off the diagonal fails %d times, "
                                      "even characters %d times" % (n1, n2))
    for ff, esp in ((7, 1), (23, 1), (31, 3), (47, 1), (11, 0), (13, 0), (29, 0), (73, 4)):
        r = semisuma(ff)
        f += _linea(r["ceros"] == cuantos(ff) == esp,
                    "f=%d prime: %d vanishing components = d/t count %d" % (ff, r["ceros"],
                                                                            cuantos(ff)))
    r21 = semisuma(21)
    f += _linea(r21["ceros"] == 2 and _declina(cuantos, 21),
                "f=21: 2 components vanish, and the count d/t declines (composite)")
    difs = [ff for ff in range(3, 400, 4) if all(k == 1 for k in factoriza(ff).values())
            and dirichlet(ff)[0] != dirichlet(ff)[1]]
    f += _linea(not difs, "Dirichlet's family, squarefree f = 3 mod 4, f < 400: no discrepancy")
    f += _linea(_declina(semisuma, 10) and _declina(semisuma, SEMISUMA_MAX + 2) and
                _declina(dirichlet, 13), "f even, f above the range, f = 1 mod 4 (Dirichlet): "
                                         "decline")

    print("the null cells of Conrey's F_q")
    c = conrey(103)
    f += _linea(c["nulas"] == [47, 51, 55] and c["igual"] and not c["prop_aplica"],
                "q=103: null cells %s = R_103; q = 7 mod 8, the proposition does not apply"
                % c["nulas"])
    c = conrey(59)
    f += _linea(c["nulas"] == [] and c["prop_aplica"], "q=59: no null cell; q = 3 mod 8, "
                                                        "the proposition applies")
    for q in (7, 23, 31, 47):
        c = conrey(q)
        f += _linea((q - 1) // 2 in c["nulas"] and c["igual"],
                    "q=%d: the central cell [(q-1)/(2q), (q+1)/(2q)] is null" % q)
    bad = [q for q in primos_hasta(2000) if q > 3 and q % 4 == 3 and not
           (lambda c: c["igual"] and c["continua"] and c["espejo"])(conrey(q))]
    f += _linea(not bad, "3 < q < 2000, q = 3 mod 4: null cells = R_q, F_q continuous, mirror")
    f += _linea(_declina(conrey, 13) and _declina(conrey, 3) and _declina(conrey, 15),
                "q = 13 (1 mod 4), q = 3, q = 15: decline")

    print()
    print("SELF-TEST: %s" % ("all correct" if f == 0 else "%d FAILURES" % f))
    return f


# ------------------------------------------------------------------------------------- the CLI

def _cli_cero(args):
    m = int(args[0])
    full = len(args) > 1 and args[1] == "full"
    s, S, how = cero(m)
    if s is None:
        print("m = %d is not admissible: %s" % (m, how))
        return 0
    neg = sorted(p for p, v in s.items() if v == -1)
    print("m = %d: a completely multiplicative eps with sum_{n<=m} n eps(n) = %d" % (m, S))
    print("   found by: %s" % how)
    print("   eps(p) = -1 at %d of the %d primes <= m; eps(p) = +1 at the others" % (len(neg), len(s)))
    if full or len(neg) <= 40:
        print("   eps(p) = -1 at: %s" % neg)
    else:
        print("   first of them: %s ...   (`cero %d full` prints all)" % (neg[:30], m))
    print("   the sum was recomputed from these signs, n by n: it is 0")
    return 0


def _cli_certificado(args):
    if args:
        r = certificado(int(args[0]))
        _eco_cert(r)
        return 1 if r["failures"] else 0
    res, nf = certificados(eco=_eco_cert)
    worst = max(res, key=lambda r: r["ratio"])
    print()
    print("blocks: %d, covering [%d, %d]; pairs: %d; failures: %d" % (
        len(res), res[0]["a"], res[-1]["b"], sum(r["pairs"] for r in res), nf))
    print("largest B*(b)^2/D^2 bound: %.4f, at block %d (it must be < 1)" % (
        float(worst["ratio"]), worst["j"]))
    print("prefixes on {2,3,5,7} (factor 16), every block with b >= 60000: %s" % all(
        r["prefix16"] for r in res if r["b"] >= 60000))
    return 1 if nf else 0


def _eco_cert(r):
    print("block %3d [%9d, %9d] |Q|=%7d pairs=%2d q0=%9d alpha=%12d D=%15d  B*^2/D^2 <= %.4f  %s"
          % (r["j"], r["a"], r["b"], r["N"], r["pairs"], r["q0"], r["alpha"], r["D"],
             float(r["ratio"]), "ok" if not r["failures"] else "*** " + "; ".join(r["failures"])))
    sys.stdout.flush()


def _cli_digitos(args):
    q = int(args[0])
    s = int(args[1]) if len(args) > 1 else 1
    X = int(float(args[2])) if len(args) > 2 else 10 ** 6
    d = digitos(q, s, X)
    print("q = %d (q = %d mod 8), s = chi(q) = %+d" % (q, q % 8, s))
    print("   T_q = %d%s" % (d["T_q"], "   h(-q) = %d" % d["h"] if d["h"] is not None else ""))
    print("   R_q = %s" % d["R"])
    z = d["ceros"]
    print("   zeros of S_chi in [1, %d]: %d%s" % (X, len(z), "" if not z else
                                                  "   %s%s" % (z[:40], " ..." if len(z) > 40
                                                               else "")))
    print("   each zero: direct summation AND the digit recursion give 0; the recursion was also "
          "checked against the direct sum at every m <= %d" % d["recursion_checked_at"])
    print("   digit tree (all digits in R_q) below the bound: %d numbers, all among the zeros"
          % len(d["arbol"]))
    for texto, ok in d["partes"]:
        print("   %s: %s" % (texto, "yes" if ok else "NO"))
    for x in d["fallos"]:
        print("   *** %s ***" % x)
    return 1 if d["fallos"] else 0


def _cli_semisuma(args):
    f = int(args[0])
    r = semisuma(f)
    print("f = %d: %d odd characters (values in Z[zeta_%d])" % (f, r["caracteres"], r["N"]))
    print("   closed form 2 chi(2) T = (sum a chi(a)) (1 - chi(2)), exactly: %d of %d"
          % (r["forma_ok"], r["caracteres"]))
    print("   primitive: %d; on them T = 0 iff chi(2) = 1 holds in %d" % (
        r["primitivos"], r["crit_prim_ok"]))
    print("   T = 0 at %d characters (%d of them with chi(2) = 1)" % (r["ceros"], r["ceros_chi2_1"]))
    if es_primo(f) and f >= 5:
        print("   f prime: the count of Cor. 'how many components vanish' is %d" % cuantos(f))
    else:
        print("   f composite: the count d/t is not claimed (imprimitive characters can add zeros)")
    if f % 4 == 3 and all(k == 1 for k in factoriza(f).values()):
        T, pred = dirichlet(f)
        print("   quadratic character (.|f): half sum %d, Dirichlet's family predicts %s" % (T, pred))
    return 0 if r["forma_ok"] == r["caracteres"] and r["crit_prim_ok"] == r["primitivos"] else 1


def _cli_conrey(args):
    q = int(args[0])
    c = conrey(q)
    print("q = %d (q = %d mod 8), h(-q) = %d" % (q, q % 8, c["h"]))
    print("   F_q(y) = T(r) + y (h - C(r)), r = floor(y): continuous %s" % c["continua"])
    if c["nulas"]:
        print("   F_q vanishes identically on the cells r = %s, i.e. f_q on x in %s" % (
            c["nulas"], ", ".join("[%s, %s]" % x for x in c["celdas_x"])))
    else:
        print("   F_q vanishes identically on no cell")
    print("   those cells are exactly R_q = %s: %s" % (c["R"], c["igual"]))
    if c["prop_aplica"]:
        print("   q = 3 mod 8: Prop. 'positivity empties the digit set' applies -- if the series of")
        print("   chi_{q,-1} is >= 0 on [0, 1/2], then R_q is empty.  Here R_q is %s."
              % ("empty" if not c["R"] else "NOT empty"))
    else:
        print("   q = 7 mod 8: the proposition's hypothesis q = 3 mod 8 is not met; the central cell")
        print("   (q-1)/2 is null by Cor. 'zeros from digits' (b)")
    return 0 if c["igual"] and c["continua"] else 1


def main(argv):
    if len(argv) < 2 or argv[1] in ("-h", "--help"):
        print(__doc__)
        return 0
    cmd, args = argv[1], argv[2:]
    tabla = {"cero": _cli_cero, "certificado": _cli_certificado, "digitos": _cli_digitos,
             "semisuma": _cli_semisuma, "conrey": _cli_conrey}
    if cmd == "autotest":
        # A crash is not a pass: it is reported as a failure, with the SELF-TEST line.
        try:
            return 1 if autotest() else 0
        except Exception as e:
            print("   *** the self-test crashed: %s: %s ***" % (type(e).__name__, e))
            print()
            print("SELF-TEST: FAILED (crash)")
            return 1
    if cmd not in tabla:
        print("unknown command: %s" % cmd)
        return 2
    try:
        return tabla[cmd](args)
    except ValueError as e:
        print("declined: %s" % e)
        return 2
    except IndexError:
        print("missing argument; see  python segmento_II.py --help")
        return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv))
