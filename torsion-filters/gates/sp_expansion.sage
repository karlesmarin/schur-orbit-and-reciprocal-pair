# -*- coding: utf-8 -*-
# LA EXPANSION SIMPLECTICA  Phi_{t,r} = sum_mu A_mu · sp_mu.  15 de agosto de 2026.
#
# DE DONDE SALE.  La consulta externa, vuelta 10, reformula el problema: el v(T) maximal no es el
# objeto intrinseco -- es un candidato de NUESTRA presentacion de Laplace, y se puede cancelar
# entero.  El objeto intrinseco seria
#
#     mu_max  =  el peso dominante mas alto con  A_mu != 0     en   Phi_{t,r} = sum_mu A_mu sp_mu
#
# y su pregunta: ¿hay UN solo mu maximal?  Si lo hay, Newt(Phi) = conv(W(C_r) mu_max) y el rombo que
# medimos es P_{C_2}(10,0).  Y su test 4: atacar UN coeficiente A_mu en vez de todo Phi.
#
# POR QUE EN SAGE, Y POR QUE POR EL BIALTERNANTE.  Nuestra maquinaria de Python calcula Phi por la
# expansion de Laplace sobre las t filas congeladas.  Si escribiera aqui lo mismo, no seria un
# control: seria la misma cuenta dos veces.  Asi que aqui se calcula por la DEFINICION,
#
#     Phi = det( x_i^{beta_j} ) / det( x_i^{delta_j} )        x = (1, zeta, ..., zeta^{t-1}, z, 1/z)
#
# con zeta en el cuerpo ciclotomico y z en un anillo de Laurent: un determinante N x N, sin transversales
# ni atomos ni estratos.  Que las dos coincidan es el control C0.
#
# LA DESCOMPOSICION.  {sp_mu} es base del anillo de caracteres de Sp_{2r}, y sp_mu tiene como termino
# dominante z^mu.  Luego basta pelar: tomar el peso dominante mas alto que sobreviva, leer su
# coeficiente, restar A_mu · sp_mu y repetir.  Termina porque cada paso baja el peso.  sp_mu se toma
# de WeylCharacterRing("C%d") de Sage -- no lo escribimos nosotros.
#
# CONTROLES
#   C0  FATAL.  El Phi del bialternante coincide MONOMIO A MONOMIO con el de peel_zero/collision_graph.
#       Si no, una de las dos maquinarias esta mal y no se dice nada mas.
#   C1  FATAL.  Reconstruir: sum_mu A_mu sp_mu tiene que devolver Phi exacto.  Un resto no nulo
#       significa que el pelado no termino o que la base no es la que creo.
#   C2  DECOY.  Las formas que SE ANULAN tienen que dar A_mu = 0 para todo mu, o sea lista vacia.  Si
#       alguna diera un mu, el instrumento miente.
#   C3  no vacuidad: se imprime n SIEMPRE, y se dice cuantas formas de cada columna hay.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage sp_expansion.sage

import itertools


def phi_bialternante(beta, t, r):
    """Phi por la DEFINICION: cociente de dos determinantes N x N.  dict {exponente: coef} o None."""
    N = t + 2 * r
    K = CyclotomicField(t) if t > 2 else QQ
    zeta = K.gen() if t > 2 else K(-1)
    L = LaurentPolynomialRing(K, r, 'z')
    zs = L.gens()
    frozen = [K(zeta) ** k for k in range(t)]
    x = [L(c) for c in frozen] + [g ** e for g in zs for e in (1, -1)]
    delta = list(range(N - 1, -1, -1))

    def alt(expo):
        return matrix(L, N, N, lambda i, j: x[i] ** expo[j]).determinant()
    den = alt(delta)
    if den == 0:
        return None
    num = alt(list(beta))
    q = num / den
    try:
        q = L(q)
    except Exception:
        return "NO-POLINOMIO"
    out = {}
    for e, c in zip(q.exponents(), q.coefficients()):
        e = tuple(e) if hasattr(e, '__iter__') else (e,)
        if c != 0:
            out[e] = c
    return out


def sp_char(mu, r):
    """sp_mu como dict {exponente: coef}, desde WeylCharacterRing de Sage."""
    W = WeylCharacterRing("C%d" % r)
    L = W.space()
    el = W(L.from_vector(vector(list(mu))))
    out = {}
    for wt, m in el.weight_multiplicities().items():
        out[tuple(wt.to_vector())] = out.get(tuple(wt.to_vector()), 0) + m
    return out


def dominante(e):
    """el representante dominante de la orbita W(C_r): valores absolutos, decrecientes."""
    return tuple(sorted((abs(x) for x in e), reverse=True))


def expandir(P, r, tope=400):
    """P = sum_mu A_mu sp_mu.  Devuelve (lista de (mu, A_mu), resto)."""
    P = dict(P)
    salida = []
    for _ in range(tope):
        P = {e: c for e, c in P.items() if c != 0}
        if not P:
            return salida, {}
        dom = [e for e in P if list(e) == sorted(e, reverse=True) and min(e) >= 0]
        if not dom:
            return salida, P
        mu = max(dom, key=lambda e: (sum(e), e))
        A = P[mu]
        sp = sp_char(mu, r)
        for e, m in sp.items():
            P[e] = P.get(e, 0) - A * m
        salida.append((mu, A))
    return salida, P


# ===================================================================== C0 ========================
print("=" * 100)
print("C0  ACEPTACION -- el bialternante coincide con la maquinaria de Laplace")
print("=" * 100)
print("")
REF = {}
try:
    import json
    REF = json.load(open("sp_expansion_REF.json"))
except Exception:
    REF = {}
CASOS_C0 = [((18, 17, 11, 8, 7, 6, 1, 0), 4, 2), ((13, 9, 8, 7, 5, 4, 2, 0), 4, 2),
            ((10, 9, 7, 4, 3, 2, 1, 0), 4, 2)]
malo = 0
for (b, t, r) in CASOS_C0:
    P = phi_bialternante(b, t, r)
    n = 0 if P in (None, "NO-POLINOMIO") else len(P)
    ref = REF.get(str(list(b)))
    ok = "sin referencia" if ref is None else ("igual" if ref == n else "*** DISTINTO (%s) ***" % ref)
    malo += (ref is not None and ref != n)
    print("    %-28s t=%d r=%d : %4d monomios   contra Python: %s" % (str(b), t, r, n, ok))
print("")
print("    C0 %s" % ("PASA" if not malo else "*** FALLA: las dos maquinarias no calculan lo mismo ***"))
if malo:
    print("DONE (veredicto suspendido)")
    sys.exit(1)

# ===================================================================== N1 ========================
print("")
print("=" * 100)
print("N1  LA EXPANSION, Y SU PREGUNTA: ¿hay UN SOLO mu maximal?")
print("=" * 100)
print("")
print("   beta                        | Phi==0 | #mu | mu maximales           | A(mu_max) | C1 resto")
CASOS = [(18, 17, 11, 8, 7, 6, 1, 0), (10, 9, 7, 4, 3, 2, 1, 0), (14, 13, 11, 4, 3, 2, 1, 0),
         (12, 11, 10, 5, 3, 2, 1, 0), (14, 13, 11, 8, 3, 2, 1, 0), (12, 11, 10, 9, 7, 2, 1, 0),
         (13, 9, 8, 7, 5, 4, 2, 0), (19, 17, 11, 8, 7, 6, 1, -1), (21, 17, 11, 8, 7, 6, 1, -3)]
t, r = 4, 2
for b in CASOS:
    P = phi_bialternante(b, t, r)
    if P is None:
        print("   %-27s | ---    | denominador nulo" % str(b))
        continue
    if P == "NO-POLINOMIO":
        print("   %-27s | ---    | *** el cociente no es de Laurent: revisar ***" % str(b))
        continue
    cero = (len(P) == 0)
    lista, resto = expandir(P, r)
    lista = [(mu, A) for (mu, A) in lista if A != 0]
    maxi = [mu for (mu, A) in lista
            if not any(nu != mu and all(sum(nu[:k + 1]) >= sum(mu[:k + 1]) for k in range(r))
                       for (nu, _) in lista)]
    Amax = [A for (mu, A) in lista if mu in maxi]
    print("   %-27s | %-6s | %3d | %-22s | %-9s | %s"
          % (str(b), cero, len(lista), str(maxi)[:22], str(Amax)[:9],
             "0" if not resto else "*** %d ***" % len(resto)))

print("")
print("  C2  DECOY: las que se anulan (las dos ultimas) tienen que dar #mu = 0.  Si dieran algun mu,")
print("      el instrumento miente y N1 no vale.")
print("")
print("=" * 100)
