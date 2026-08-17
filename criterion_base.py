# -*- coding: utf-8 -*-
"""Las rutinas compartidas por los cinco guiones del criterio de la Seccion 8.

POR QUE EXISTE ESTE FICHERO.  `criterion_S.py`, `criterion_control.py`,
`criterion_sage_check.sage`, `sign_lemma.py` y `step_law.py` importaban estas nueve funciones de
dos guiones mucho mayores que pertenecen a una linea de trabajo distinta, todavia sin publicar.
Distribuir aquellos habria supuesto publicar de rebote esa linea; no distribuir nada habria dejado
cinco resultados citados en las Secciones 8 y 9 sin codigo detras.  Asi que aqui estan las
funciones que los cinco usan, y solo esas.

Los cuerpos son los originales, extraidos automaticamente y verificados caracter a caracter: este
fichero no es una reescritura, es un recorte.

  perm_sign, setup, deg_of, all_transversals   la expansion de Laplace sobre las transversales
  split_sign, alt, stratify, stratum, measure  la filtracion por grados y el primer estrato no nulo

Autores: Carles Marin + Claude (AI assistant).
"""
import itertools
from collections import defaultdict
from itertools import combinations, permutations


# --- de second_stratum.py ---
def perm_sign(seq):
    """signo de la permutacion que ordena seq (elementos distintos)."""
    s = 1
    a = list(seq)
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            if a[i] > a[j]:
                s = -s
    return s

# --- de second_stratum.py ---
def setup(beta, t):
    """clases de residuos; None si alguna esta vacia (hipotesis de ocupacion (O))."""
    cl = defaultdict(list)
    for i, b in enumerate(beta):
        cl[b % t].append(i)
    if len(cl) < t:
        return None
    E = sorted(k for k in cl if len(cl[k]) >= 2)
    C = {k: sorted((beta[i] for i in cl[k]), reverse=True) for k in E}
    return dict(cl), E, C

# --- de second_stratum.py ---
def deg_of(T, r):
    return sum(T[:r]) - sum(T[r:])

# --- de second_stratum.py ---
def all_transversals(beta, cl, r, t):
    """[(clave por clase, T, w, deg)] sobre todas las transversales."""
    N = len(beta)
    out = []
    keys = sorted(cl)
    for pick in itertools.product(*[cl[k] for k in keys]):
        P = sorted(pick)
        Ps = set(P)
        T = tuple(beta[i] for i in range(N) if i not in Ps)
        w = perm_sign([beta[i] % t for i in P])
        if sum(P) % 2:
            w = -w
        sel = {k: beta[i] for k, i in zip(keys, pick)}
        out.append((sel, T, w, deg_of(T, r)))
    return out

# --- de depth_histogram.py ---
def split_sign(S, Sc, r):
    """signo de la permutacion que manda S a las columnas pares y S^c a las impares."""
    q = [0] * (2 * r)
    for j in range(r):
        q[S[j]] = 2 * j
        q[Sc[j]] = 2 * j + 1
    return perm_sign(q)

# --- de depth_histogram.py ---
def alt(vals, r):
    """el alternante: r! monomios con signo, clave = exponentes por variable."""
    return {tuple(vals[p[j]] for j in range(r)): perm_sign(p)
            for p in permutations(range(r))}

# --- de depth_histogram.py ---
def stratify(tm, r):
    """{grado total: [(w, T, S, S^c)]} -- los C(2r,r) bloques de Laplace de cada transversal."""
    buckets = defaultdict(list)
    idx = tuple(range(2 * r))
    for w, T in tm:
        tot = sum(T)
        for S in combinations(idx, r):
            Sc = tuple(a for a in idx if a not in S)
            buckets[2 * sum(T[a] for a in S) - tot].append((w, T, S, Sc))
    return buckets

# --- de depth_histogram.py ---
def stratum(bucket, r):
    """el estrato entero: {monomio: coeficiente}, ya sin los ceros."""
    acc = defaultdict(int)
    for (w, T, S, Sc) in bucket:
        A = alt([T[a] for a in S], r)
        B = alt([-T[a] for a in Sc], r)
        base = w * split_sign(S, Sc, r)
        for ka, ca in A.items():
            c = base * ca
            for kb, cb in B.items():
                acc[tuple(ka[j] + kb[j] for j in range(r))] += c * cb
    return {k: v for k, v in acc.items() if v}

# --- de depth_histogram.py ---
def measure(tm, r, want=1):
    """None si Phi_t == 0; si no (primer grado no nulo, grados que SE CANCELAN por encima,
    grados no nulos hasta 'want', conjunto de grados con soporte).

    CONVENIO DE ESTRATO -- el de depth.py:205 y el de depth6_check, que es el de los enunciados:
    los estratos son D1, D1-2, D1-4, ... se cuenten o no.  Un grado por encima del primero no nulo
    puede estar vacio por DOS razones distintas, y aqui se separan porque no son lo mismo para
    quien quiera probar una cota: o tiene monomios que se cancelan, o no tiene monomio ninguno."""
    B = stratify(tm, r)
    first = None
    spectrum = []
    cancel = []
    for s in sorted(B, reverse=True):
        if stratum(B[s], r):
            if first is None:
                first = s
            spectrum.append(s)
            if len(spectrum) >= want:
                break
        elif first is None:
            cancel.append(s)
    if first is None:
        return None
    return first, cancel, spectrum, set(B)
