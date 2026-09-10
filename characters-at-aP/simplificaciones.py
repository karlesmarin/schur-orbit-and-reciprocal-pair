# -*- coding: utf-8 -*-
# SIMPLIFICAR LAS FORMULAS, Y VER QUE QUEDA.   8 de septiembre de 2026.
#
# Con el criterio, el valor con signo, la factorizacion y los recuentos ya en la mano, hay tres
# simplificaciones que se ven a ojo.  Esta gate las mide.
#
# S1  EL RECUENTO, EN UNA SOLA FRASE.  En regimen regular medimos
#         q impar -> prod (q - (2i-1))      q par -> prod (q - 2i).
#     Pero en C_m los EXPONENTES son 1,3,...,2m-1 y los GRADOS son 2,4,...,2m.  O sea:
#         q impar -> prod (q - m_i)   [exponentes]      q par -> prod (q - d_i)   [grados]
#     Si eso vale, la paridad de q no es un caso feo: es lo que decide si el sistema se lee por
#     exponentes o por grados.  Se mide en C, y de paso en B, D, G2, F4, E6, para ver hasta donde.
#
# S2  EL PERFIL DE OCUPACION, CUANDO d ES IMPAR.  Con d impar las capacidades se saturan, luego el
#     perfil de ocupacion de CUALQUIER superviviente deberia ser el de rho.  Si es asi, el grupo
#     producto de la factorizacion es el MISMO para todo lambda superviviente --- lambda solo elige
#     el peso maximo ---, y ese grupo tiene el perfil de bloques del centralizador del elemento.
#     Eso convierte la factorizacion en una descomposicion indexada por bloques, no en una
#     coincidencia numerica.
#
# S3  EL SIGNO, ¿SE SIMPLIFICA?  E = sum_a (floor(p(ell,a)/q) - floor(p(rho,a)/q)).  Escribiendo
#     (x,a) = q*floor + resto y usando sum_{a>0}(x,a) = 2(x,rho), sale
#         E = [2(lambda,rho)*p - sum_a(resto_ell) + sum_a(resto_rho)] / q.
#     Y (B) dice que los restos coinciden SALVO r <-> -r.  Cada pareja volteada cambia la suma en
#     q - 2r.  Asi que E deberia ser  2p(lambda,rho)/q  corregido por los volteos.  Se mide contra
#     el candidato limpio y contra el corregido.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python simplificaciones.py > simplificaciones_OUT.txt 2>&1

import io
import json
import sys
from itertools import product as iproduct
from math import gcd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SEP = "=" * 96


def prod(v):
    r = 1
    for x in v:
        r *= x
    return r


def dominants(m, top):
    def rec(k, hi):
        if k == 0:
            yield ()
            return
        for v in range(hi, -1, -1):
            for rest in rec(k - 1, v):
                yield (v,) + rest
    return rec(m, top)


def pares_C(v):
    m = len(v)
    out = []
    for i in range(m):
        for j in range(i + 1, m):
            out.append(v[i] - v[j]); out.append(v[i] + v[j])
    for i in range(m):
        out.append(2 * v[i])
    return out


def cls(x, q):
    return min(x % q, (-x) % q)


# ---------------------------------------------------------------- S1
print(SEP)
print("S1 -- ¿EXPONENTES SI q ES IMPAR, GRADOS SI ES PAR?   (grado = exponente + 1)")
print(SEP)
DAT = json.load(open("root_data.json", encoding="utf-8"))
print(f"{'tipo':>5} {'q':>4} {'S(q)':>10} {'prod(q-m_i)':>12} {'prod(q-d_i)':>12} {'cual':>10}")
res = {"exp_impar": 0, "gra_par": 0, "mal_impar": 0, "mal_par": 0}
for d in DAT:
    n, R, exps = d["rango"], d["raices"], d["exponentes"]
    if n > 4:
        continue
    grados = [m + 1 for m in exps]
    hmax = max(sum(v) * c for c, v in R)
    for q in range(2 * n + 1, 2 * n + 9):     # regimen regular
        if (q ** n) * len(R) > 30_000_000:
            continue
        r_q = sum(1 for c, v in R if (c * sum(v)) % q == 0)
        S = 0
        for x in iproduct(range(q), repeat=n):
            k = 0
            for c, v in R:
                s = 0
                for j in range(n):
                    s += x[j] * v[j]
                if (c * s) % q == 0:
                    k += 1
                    if k > r_q:
                        break
            if k == r_q:
                S += 1
        pe = prod([q - m for m in exps])
        pg = prod([q - g for g in grados])
        cual = ("exponentes" if S == pe else ("grados" if S == pg else "ninguna"))
        if q % 2:
            res["exp_impar" if S == pe else "mal_impar"] += 1
        else:
            res["gra_par" if S == pg else "mal_par"] += 1
        print(f"{d['tipo']:>5} {q:>4} {S:>10} {pe:>12} {pg:>12} {cual:>10}")
print()
print(f"   q impar y sale por EXPONENTES: {res['exp_impar']}   no: {res['mal_impar']}")
print(f"   q par   y sale por GRADOS ...: {res['gra_par']}   no: {res['mal_par']}")

# ---------------------------------------------------------------- S2
print(SEP)
print("S2 -- CON d IMPAR, ¿todo superviviente tiene el perfil de ocupacion de rho?")
tot = mal = 0
ejemplo = None
for m in range(2, 10):
    t = 2 * m + 2
    rho = [m - i for i in range(m)]
    Prho = pares_C(rho)
    top = 9 if m <= 5 else 6
    for k in range(1, t):
        dd = gcd(k, t); q = t // dd
        if q < 2 or dd % 2 == 0:
            continue
        r_q = sum(1 for b in Prho if b % q == 0)
        perfil_rho = tuple(sorted((c, sum(1 for x in rho if cls(x, q) == c))
                                  for c in range(q // 2 + 1)))
        for lam in dominants(m, top):
            ell = [lam[i] + rho[i] for i in range(m)]
            if sum(1 for a in pares_C(ell) if a % q == 0) != r_q:
                continue
            perfil = tuple(sorted((c, sum(1 for x in ell if cls(x, q) == c))
                                  for c in range(q // 2 + 1)))
            tot += 1
            if perfil != perfil_rho:
                mal += 1
                if ejemplo is None:
                    ejemplo = (m, k, q, lam, perfil, perfil_rho)
print(f"   supervivientes con d impar: {tot};  con perfil distinto del de rho: {mal}")
if ejemplo:
    print("   ejemplo:", ejemplo)
print("   " + ("=> el perfil es el de rho SIEMPRE: el grupo de la factorizacion no depende de lambda"
               if mal == 0 else "=> el perfil NO es siempre el de rho"))

# ---------------------------------------------------------------- S3
print(SEP)
print("S3 -- ¿SE SIMPLIFICA EL SIGNO?   candidato limpio:  E == 2p(lambda,rho)/q   (mod 2)")
print(f"{'m':>3} {'k':>3} {'d':>3} {'q':>4} {'vivos':>7} {'E==cand':>9} {'E-cand par':>11}")
T = {"ok": 0, "bad": 0}
difs = {}
for m in range(2, 8):
    t = 2 * m + 2
    rho = [m - i for i in range(m)]
    Prho = pares_C(rho)
    top = 8 if m <= 4 else 6
    for k in range(1, t):
        dd = gcd(k, t); q = t // dd; p = k // dd
        if q < 2:
            continue
        r_q = sum(1 for b in Prho if b % q == 0)
        Erho = sum((p * b) // q for b in Prho)
        vivos = ok = 0
        for lam in dominants(m, top):
            ell = [lam[i] + rho[i] for i in range(m)]
            P = pares_C(ell)
            if sum(1 for a in P if a % q == 0) != r_q:
                continue
            vivos += 1
            E = sum((p * a) // q for a in P) - Erho
            lr = sum(lam[i] * rho[i] for i in range(m))     # (lambda, rho)
            cand = (2 * p * lr) // q if (2 * p * lr) % q == 0 else None
            if cand is not None and (E - cand) % 2 == 0:
                ok += 1
            difs[(E - cand) % 2 if cand is not None else "no-entero"] = \
                difs.get((E - cand) % 2 if cand is not None else "no-entero", 0) + 1
        if vivos == 0:
            continue
        T["ok"] += ok; T["bad"] += vivos - ok
        print(f"{m:>3} {k:>3} {dd:>3} {q:>4} {vivos:>7} {ok:>5}/{vivos:<4}"
              + ("   <-- no cuadra" if ok != vivos else ""))
print()
print(f"   E == 2p(lambda,rho)/q (mod 2):  {T['ok']} si, {T['bad']} no")
print("   reparto de (E - candidato) mod 2:", difs)
print(SEP)
