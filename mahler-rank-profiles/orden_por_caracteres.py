# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marin. All rights reserved.
# Carles Marin <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
#
# QUE: calculo del orden (piso 1, a = M = 1) por caracteres, sin construir el anillo (1-oct-2026).
#      Caso semisimple p no divide n phi(n). En T (x) F_p-barra las componentes isotipicas se multiplican caracter a
#      caracter (L_chi L_psi = L_{chi psi}), asi que gr(A) queda descrito por SOPORTES S_d (caracteres presentes en el
#      grado d):
#        d = 0: {trivial};  1 <= d <= p-2: vacio;  d = p-1: {trivial};
#        d = p+j (0 <= j <= p-2), ventana exacta (Teorema A, CRITERIO_N3_PRUEBA.md §1 + Lema de Abel):
#            chi no trivial de paridad (-1)^{j+1}: presente sii (1 - 2^j chi(2)) B_{j+1, chi_prim} != 0 mod P;
#            trivial (j impar): presente sii (1 - 2^j) B_{j+1} != 0 mod p, o j = p-2;
#        d = 2p-1: soporte de V_0 (Teorema 1);
#        d >= 2p: cierre por productos, S_{a+b} contiene {chi psi : chi en S_a, psi en S_b} (punto fijo).
#      El cierre da un SUBANILLO de gr(A): cota inferior de las capas, luego cota superior de e (certificado [E, 2E)).
#      Las capas vacias en [p-1, 2p-1] (exactas) dan la cota inferior. Si coinciden, e exacto.
#      Validacion: e medido por el calculo directo (conjetura_C_{lote,p17}_OUT.txt).
# FECHA: 1-oct-2026
# Uso (docker, /work): sage -python orden_por_caracteres.py
import re

from sage.all import GF, QQ, DirichletGroup, Integers, bernoulli, bernoulli_polynomial, euler_phi, lcm


def B_prim(chi, k, F):
    f = chi.conductor()
    if f == 1:
        return F(QQ(bernoulli(k)))
    Bk = bernoulli_polynomial(QQ(0), k)
    s = F(0)
    for a in range(1, f + 1):
        ca = chi(a)
        if ca:
            s += ca * F(QQ(f) ** (k - 1) * (bernoulli_polynomial(QQ(a) / f, k) - Bk))
    return s


def orden(n, p, dmax_factor=4):
    lam = Integers(n).unit_group_exponent()
    f = 1
    while (p ** f - 1) % lam:
        f += 1
    F = GF(p ** f, 'g')
    z = F.multiplicative_generator() ** ((p ** f - 1) // lam)
    G = list(DirichletGroup(n, F, zeta=z, zeta_order=lam))
    idx = {chi: i for i, chi in enumerate(G)}
    triv = next(i for i, c in enumerate(G) if c.is_trivial())
    mult = [[idx[G[i] * G[j]] for j in range(len(G))] for i in range(len(G))]
    impares = {i for i, c in enumerate(G) if c(-1) == -1}
    pares = set(range(len(G))) - impares
    d = euler_phi(n) // 2
    Dmax = dmax_factor * p + 4
    S = [set() for _ in range(Dmax + 1)]
    S[0] = {triv}
    S[p - 1] = {triv}
    prim = {i: G[i].primitive_character() for i in range(len(G))}
    for j in range(p - 1):
        dd = p + j
        par = (-1) ** (j + 1)
        for i, c in enumerate(G):
            if c(-1) != par:
                continue
            if i == triv:
                if j % 2 == 1 and (j == p - 2 or F((1 - 2 ** j) * QQ(bernoulli(j + 1))) != 0):
                    S[dd].add(i)
                continue
            if (1 - F(2) ** j * c(2)) * B_prim(prim[i], j + 1, F) != 0:
                S[dd].add(i)
    S[2 * p - 1] = set(S[p])  # Teorema 1: G_{2p-1} = V_0
    # cierre por productos (punto fijo)
    cambio = True
    while cambio:
        cambio = False
        for a in range(p - 1, Dmax + 1):
            if not S[a]:
                continue
            for b in range(a, Dmax + 1 - a):
                if not S[b]:
                    continue
                nuevo = {mult[x][y] for x in S[a] for y in S[b]}
                if not nuevo <= S[a + b]:
                    S[a + b] |= nuevo
                    cambio = True
    lleno = [len(S[k]) == (len(pares) if k % 2 == 0 else len(impares)) for k in range(Dmax + 1)]
    # cota superior: primer E con [E, 2E) lleno
    E_sup = next((E for E in range(1, Dmax // 2 + 1) if all(lleno[k] for k in range(E, 2 * E))), None)
    # cota inferior: ultima capa no llena en la zona exacta [p-1, 2p-1], +1 ; y como minimo p (d > 1) o p-1 (d = 1)
    no_llenas = [k for k in range(p - 1, 2 * p) if not lleno[k]]
    E_inf = max(no_llenas) + 1 if no_llenas else (p - 1 if d == 1 else p)
    return E_inf, E_sup, d


PAT = re.compile(r"n=\s*(\d+) p=\s*(\d+) d=\s*(\d+) .*?e=(\d+|None)")
casos = []
for nom in ('conjetura_C_lote_OUT.txt', 'conjetura_C_p17_OUT.txt'):
    for ln in open('/work/' + nom, encoding='utf-8'):
        m = PAT.search(ln)
        if m and m[4] != 'None':
            casos.append((int(m[1]), int(m[2]), int(m[4])))
tot = exacto = dentro = 0
malos = []
for n, p, e in casos:
    if euler_phi(n) % p == 0:
        continue
    lo, hi, d = orden(n, p)
    tot += 1
    if lo == hi == e:
        exacto += 1
    if hi is not None and lo <= e <= hi:
        dentro += 1
    else:
        malos.append((n, p, e, lo, hi))
    print("n=%2d p=%2d d=%2d  e medido %2d   herramienta [%s, %s]  %s" % (
        n, p, d, e, lo, hi, 'EXACTO' if lo == hi == e else ('dentro' if hi and lo <= e <= hi else 'FUERA')),
        flush=True)
print("casos %d: e exacto %d, e dentro del intervalo %d; fuera: %s" % (tot, exacto, dentro, malos))
