# -*- coding: utf-8 -*-
# LA DIVISION POR Delta_t COMO OPERADOR DE DIFERENCIAS DIVIDIDAS.   16 de agosto de 2026.
#
# Lo que queda de (L1) es una sola frase: la sustitucion hacia atras por
#
#       Delta_t = prod_j ( z_j^{t/2} - z_j^{-t/2} )
#
# no pasa nunca de 1.  Hasta ahora eso se comprobaba PELANDO (algoritmo, no formula).  Aqui se
# escribe la division como una FORMULA CERRADA, siguiendo la ruta clasica:
# operadores de diferencias divididas (Demazure), en la forma que les da Harada-Landweber-Sjamaar,
# arXiv:0906.1629 -- alli  d_alpha(u) = (u - e^{-alpha} s_alpha u)/(1 - e^{-alpha})  y la composicion
# sobre w_0 es  J(u)/d :  antisimetrizar y dividir por el denominador de Weyl.
#
# LA OBSERVACION.  Delta_t ES el denominador de Weyl de un sistema de raices: A_1^r con la raiz
# escalada por t (raices  alpha_j = t f_j ).  Su grupo de Weyl son los r cambios de signo, que
# CONMUTAN, asi que su  d_{w_0}  es el producto de r operadores de rango uno.  Y un operador de
# Demazure de rango uno es una SUMA A LO LARGO DE UNA PROGRESION:
#
#       pi(z^a) = z^a + z^{a-alpha} + ... + z^{-a}.
#
# En las coordenadas dobladas  X = 2(mu + rho_{D_r})  -- las mismas en que multiplicar por Delta_t
# es  X -> X + t*eps  con signo prod(eps) -- la inversa se escribe entera:
#
#       c(X)  =  sum over k in (2Z_{>=0}+1)^r  of  nu~( X + t*k ),
#
# con nu~ la extension de nu a todo Z^r por antisimetria de W(D_r) (signo al enderezar, 0 en la
# pared).  Rango 1: c(X) = nu(X+t) + nu(X+3t) + ... y telescopa.  La suma es finita porque nu tiene
# soporte finito, y c se anula por debajo del soporte porque la suma TOTAL sobre cada progresion es
# cero -- que es justo la divisibilidad.
#
# POR QUE IMPORTA.  nu ya es {0,+-1} con soporte descrito por transversales (prop:transversal).  Con
# la formula, (L1) deja de ser "una division nunca produce un 2" y pasa a ser:
#
#       la suma con signo de nu a lo largo de una progresion de paso 2t esta en {0,+-1}.
#
# Puro conteo.  Y hay un mecanismo candidato inmediato: que a lo sumo UN termino de cada progresion
# sea no nulo.  Eso lo mide C3, y si sale, (L1) queda probada modulo un lema de disjuncion.
#
# CONTROLES
#   C0 (FATAL)  la formula cerrada == la c del pelado, peso a peso.
#   C1          vuelta atras: multiplicar la c de la formula por Delta_t devuelve nu.
#   C2          max |c| sobre todo.
#   C3          histograma de cuantos terminos NO NULOS tiene cada progresion.  <-- el contenido
#   S1 senuelo  progresion par (2j) en vez de impar (2j+1).
#   S2 senuelo  sin el signo del enderezado.
#   S3 senuelo  solo el termino k=(1,...,1) -- la lectura ingenua de una esquina.
#
# Los ayudantes plegar/sgn_perm/eps_t/delta_dec/enderezar_D/desplazar/nu_de vienen de
# _probe_quiralidad.py, donde quedaron validados (Q1 antisimetria, Q3 654/654).
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python divided_differences.py        (aritmetica entera, no necesita Sage)

import itertools
import json
from collections import Counter, defaultdict

CASOS = [(3, 2, 6), (5, 2, 4), (7, 2, 3), (3, 3, 4)]


def plegar(v, t):
    v %= t
    if v == 0:
        return (0, 1)
    if 2 * v == t:
        return (t // 2, 1)
    return (v, 1) if v < t - v else (t - v, -1)


def sgn_perm(perm):
    n, s, visto = len(perm), 1, [False] * len(perm)
    for i in range(n):
        if visto[i]:
            continue
        j, L = i, 0
        while not visto[j]:
            visto[j] = True
            j = perm[j]
            L += 1
        if L % 2 == 0:
            s = -s
    return s


def jacobi(a, n):
    a %= n
    r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            r = -r
        a %= n
    return r if n == 1 else 0


def eps_t(t, mp):
    return jacobi((-2) % t, t) ** ((t + 3) // 2) * (1 if (mp * (mp - 1) // 2) % 2 == 0 else -1)


def delta_dec(A, t, mp):
    cl, ep = [], []
    for v in A:
        c_, e_ = plegar(v, t)
        cl.append(c_)
        ep.append(e_)
    if sorted(cl) != list(range(1, mp + 1)):
        return 0
    s = sgn_perm([mp - c for c in cl])
    for e in ep:
        s *= e
    return int(s)


def enderezar_D(x):
    """Forma canonica de un indice alternante de D_r.  Devuelve (canonico, signo) o None si pared."""
    r = len(x)
    a = [abs(int(v)) for v in x]
    if len(set(a)) != r:
        return None
    idx = sorted(range(r), key=lambda i: -a[i])
    s = sgn_perm(list(idx))
    y = [int(x[i]) for i in idx]
    neg = sum(1 for v in y if v < 0)
    cero = any(v == 0 for v in y)
    y = [abs(v) for v in y]
    if neg % 2 == 1 and not cero:
        y[-1] = -y[-1]
    return (tuple(y), s)


def desplazar(x, paso, r):
    """Multiplicar el caracter de indice x por prod_j (z^{paso/2} - z^{-paso/2})."""
    out = defaultdict(lambda: 0)
    for eps in itertools.product((1, -1), repeat=r):
        sg = 1
        for e in eps:
            sg *= e
        e2 = enderezar_D(tuple(int(x[j]) + paso * eps[j] for j in range(r)))
        if e2 is None:
            continue
        out[e2[0]] += sg * e2[1]
    return {k: v for k, v in out.items() if v != 0}


def cabeza(d):
    return max(d, key=lambda k: (sum(k), k))


def nu_de(Lam, t, r):
    """El numerador GKRS, en coordenadas dobladas.  prop:transversal."""
    mp, Rp = (t - 1) // 2, (t - 1) // 2 + r
    V = [2 * Lam[i] + 2 * (Rp - i) - 1 for i in range(Rp)]
    d = {}
    for i, x in enumerate(V):
        d.setdefault(plegar(x, t)[0], []).append(i)
    if any(not d.get(j) for j in range(1, mp + 1)):
        return {}
    E = eps_t(t, mp)
    out = {}
    for pick in itertools.product(*[d[j] for j in range(1, mp + 1)]):
        S = frozenset(pick)
        if len(S) != mp:
            continue
        Sc = [i for i in range(Rp) if i not in S]
        A = sorted([V[i] for i in S], reverse=True)
        dv = delta_dec(A, t, mp)
        if not dv:
            continue
        for qui in (1, -1):
            libre = sorted([V[i] for i in Sc], reverse=True)
            libre[-1] *= qui
            orden = sorted(S, key=lambda i: -V[i]) + sorted(Sc, key=lambda i: -V[i])
            sg = sgn_perm([orden.index(i) for i in range(Rp)])
            if qui == -1:
                sg = -sg
            out[tuple(libre)] = out.get(tuple(libre), 0) + sg * E * dv
    return {k: v for k, v in out.items() if v != 0}


def dividir(nu, t, r, tope=20000):
    """La REFERENCIA: division por pelado, el algoritmo que ya usaba el paper."""
    P = dict(nu)
    c = {}
    for _ in range(tope):
        P = {k: v for k, v in P.items() if v != 0}
        if not P:
            return c, {}
        y = cabeza(P)
        cand = None
        for eps in itertools.product((1, -1), repeat=r):
            e = enderezar_D(tuple(int(y[j]) - t * eps[j] for j in range(r)))
            if e is None:
                continue
            D = desplazar(e[0], t, r)
            if D and cabeza(D) == y:
                cand = (e[0], D)
                break
        if cand is None:
            return c, P
        x, D = cand
        if P[y] % D[y] != 0:
            return c, P
        cv = P[y] // D[y]
        c[x] = c.get(x, 0) + cv
        for k, v in D.items():
            nv = P.get(k, 0) - cv * v
            if nv == 0:
                P.pop(k, None)
            else:
                P[k] = nv
    return c, P


# ----------------------------------------------------------------------------------------------
# La formula cerrada
# ----------------------------------------------------------------------------------------------

def nu_extendida(nu, Y):
    """nu~ : la extension de nu a todo Z^r por antisimetria de W(D_r).  0 en las paredes."""
    e = enderezar_D(Y)
    if e is None:
        return 0
    return e[1] * nu.get(e[0], 0)


def c_formula(nu, t, r, paridad="impar", con_signo=True, solo_primero=False):
    """c(X) = sum_{k impar >= 1} nu~(X + t k).   Los tres senuelos entran por los flags."""
    if not nu:
        return {}
    M = max(max(abs(v) for v in k) for k in nu)          # cota del soporte
    out = {}
    rango = range(-M, M + 1)
    for X in itertools.product(rango, repeat=r):
        # solo indices dominantes regulares: X_1 > ... > X_{r-1} > |X_r|, que es la forma que
        # devuelve enderezar_D.  El ultimo puede ser negativo: es la quiralidad.
        if any(X[j] <= X[j + 1] for j in range(r - 2)):
            continue
        if not (X[r - 2] > abs(X[r - 1])):
            continue
        # los k por coordenada: impares (o pares, si el senuelo) hasta salir de la ventana
        ks = []
        for j in range(r):
            arranque = 1 if paridad == "impar" else 2
            lista = []
            k = arranque
            while X[j] + t * k <= M:
                if X[j] + t * k >= -M:
                    lista.append(k)
                k += 2
            if solo_primero:
                lista = lista[:1]
            ks.append(lista)
        if any(not L for L in ks):
            continue
        s = 0
        nz = 0
        cual = None
        for k in itertools.product(*ks):
            Y = tuple(X[j] + t * k[j] for j in range(r))
            v = nu_extendida(nu, Y) if con_signo else abs(nu_extendida(nu, Y))
            if v:
                nz += 1
                cual = k
                s += v
        # Se guarda TODO X con al menos un termino no nulo, incluidos los que suman 0.  Si el
        # histograma se condicionara a c != 0, una progresion de dos terminos que se cancelan seria
        # invisible -- y es justo el caso que decidiria si el lema de disjuncion es falso.
        if nz:
            out[X] = (s, nz, cual)
    return out


def solo_valores(d):
    return {k: v[0] for k, v in d.items() if v[0] != 0}


if __name__ == "__main__":
    print("=" * 104)
    print("LA DIVISION POR Delta_t COMO OPERADOR DE DIFERENCIAS DIVIDIDAS")
    print("=" * 104)
    print("")

    C0 = C0n = 0
    C1 = C1n = 0
    maxc = 0
    hist = Counter()
    hist_vivo = Counter()
    cancelan = Counter()
    kdist = Counter()
    S1 = S1n = S2 = S2n = S3 = S3n = 0
    testigos = []
    multiples = []
    por_caso = []
    descartados = 0

    for (t, r, cota) in CASOS:
        Rp = (t - 1) // 2 + r
        c0 = c0n = c1 = c1n = s1 = s1n = s2 = s2n = s3 = s3n = 0
        mx = desc = 0
        for Lam in itertools.product(range(cota + 1), repeat=Rp):
            if any(Lam[i] < Lam[i + 1] for i in range(Rp - 1)):
                continue
            nu = nu_de(list(Lam), t, r)
            if not nu:
                continue
            ref, resto = dividir(nu, t, r)
            if resto:
                desc += 1                     # el pelado no cerro: no hay referencia con que comparar,
                continue                      # y se CUENTA, para que no desaparezca en silencio
            ref = {k: v for k, v in ref.items() if v != 0}

            cf = c_formula(nu, t, r)
            c0n += 1
            if solo_valores(cf) == ref:
                c0 += 1
            elif len(testigos) < 3:
                testigos.append({"t": t, "r": r, "Lambda": list(Lam),
                                 "formula": {str(k): v[0] for k, v in cf.items()},
                                 "pelado": {str(k): v for k, v in ref.items()}})

            # C1  vuelta atras
            atras = defaultdict(lambda: 0)
            for X, (v, _, _) in cf.items():
                if v == 0:
                    continue
                for k, w in desplazar(X, t, r).items():
                    atras[k] += v * w
            atras = {k: v for k, v in atras.items() if v != 0}
            c1n += 1
            if atras == nu:
                c1 += 1

            for X, (v, nz, cual) in cf.items():
                mx = max(mx, abs(v))
                hist[nz] += 1
                if v:
                    hist_vivo[nz] += 1
                    kdist[cual] += 1
                else:
                    cancelan[nz] += 1
                if nz >= 2 and len(multiples) < 6:
                    multiples.append({"t": t, "r": r, "Lambda": list(Lam), "X": list(X),
                                      "c": v, "terminos_no_nulos": nz})

            # senuelos
            s1n += 1
            if solo_valores(c_formula(nu, t, r, paridad="par")) == ref:
                s1 += 1
            s2n += 1
            if solo_valores(c_formula(nu, t, r, con_signo=False)) == ref:
                s2 += 1
            s3n += 1
            if solo_valores(c_formula(nu, t, r, solo_primero=True)) == ref:
                s3 += 1

        C0 += c0; C0n += c0n; C1 += c1; C1n += c1n
        S1 += s1; S1n += s1n; S2 += s2; S2n += s2n; S3 += s3; S3n += s3n
        maxc = max(maxc, mx)
        descartados += desc
        print("  t=%d r=%d  Lambda_i <= %d :  formula == pelado %4d de %4d | vuelta atras %4d de %4d | "
              "max |c| = %d | pelado sin cerrar: %d"
              % (t, r, cota, c0, c0n, c1, c1n, mx, desc))
        por_caso.append({"t": t, "r": r, "cota": cota, "C0": c0, "n": c0n, "C1": c1, "max_c": mx,
                         "descartados": desc})

    print("")
    print("-" * 104)
    print("  C0  FATAL  la formula cerrada == el pelado      : %d de %d" % (C0, C0n))
    print("  C1         multiplicar de vuelta devuelve nu    : %d de %d" % (C1, C1n))
    print("  C2         max |c| sobre todo                   : %d" % maxc)
    print("")
    print("  C3  cuantos terminos NO NULOS tiene cada progresion, sobre TODO X con algun termino:")
    tot = sum(hist.values())
    for k in sorted(hist):
        print("        %2d termino(s) : %6d  (%5.2f%%)" % (k, hist[k], 100.0 * hist[k] / tot if tot else 0))
    print("        total X con algun termino no nulo : %d" % tot)
    print("")
    print("  C4  .hay progresiones de 2+ terminos que se CANCELAN a cero?  (el agujero de C3)")
    if cancelan:
        for k in sorted(cancelan):
            print("        %2d termino(s) sumando 0 : %6d" % (k, cancelan[k]))
    else:
        print("        NINGUNA.  Todo X con algun termino no nulo tiene c(X) != 0.")
    print("        X con c != 0 : %d   |   X con c == 0 y algun termino : %d"
          % (sum(hist_vivo.values()), sum(cancelan.values())))
    print("")
    print("  C5  .que k es el termino que sobrevive?  (si fuera siempre (1,..,1) la lectura ingenua valdria)")
    for k, n in sorted(kdist.items(), key=lambda kv: -kv[1])[:8]:
        print("        k=%-16s : %6d" % (str(k), n))
    print("        k distintos : %d" % len(kdist))
    print("")
    print("  SENUELOS (tienen que fallar)")
    print("  S1  progresion PAR    : %d de %d" % (S1, S1n))
    print("  S2  sin el signo del enderezado : %d de %d" % (S2, S2n))
    print("  S3  solo el termino k=(1,...,1) : %d de %d" % (S3, S3n))
    print("  pesos con el pelado sin cerrar (sin referencia) : %d" % descartados)
    if testigos:
        print("")
        print("  !! testigos de desacuerdo C0:")
        print("  " + json.dumps(testigos[0])[:600])
    if multiples:
        print("")
        print("  progresiones con 2 o mas terminos no nulos, primeros testigos:")
        for m in multiples[:4]:
            print("    t=%d r=%d Lambda=%s X=%s -> c=%+d con %d terminos"
                  % (m["t"], m["r"], m["Lambda"], m["X"], m["c"], m["terminos_no_nulos"]))
    print("")
    print("  LECTURA: si C0 y C1 salen enteros, la division tiene formula cerrada y (L1) es un recuento")
    print("  con signo sobre una progresion de paso 2t.  Si ademas C3 dice 1 termino SIEMPRE, entonces")
    print("  (L1) se sigue de nu en {0,+-1} y de un lema de disjuncion -- es decir, queda PROBADA.")

    json.dump({"por_caso": por_caso, "C0": [C0, C0n], "C1": [C1, C1n], "max_c": maxc,
               "hist_terminos": {str(k): v for k, v in sorted(hist.items())},
               "senuelos": {"par": [S1, S1n], "sin_signo": [S2, S2n], "solo_primero": [S3, S3n]},
               "descartados": descartados, "testigos": testigos, "multiples": multiples},
              open("divided_differences_DUMP.json", "w"), indent=1)
    print("")
    print("=" * 104)
    print("DONE")
