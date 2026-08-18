# -*- coding: utf-8 -*-
# ¿EN QUE BASE VIVE EL OBJETO IMPAR?   15 de agosto de 2026.
#
# EL PROBLEMA, Y LO ENCONTRE ESCRIBIENDO EL PAPER, NO MIDIENDO.  Para t PAR el bloque congelado es
#
#     mu_t = {1,-1} u {xi^{+-j}}      ->   alfabeto = coset de reflexion de O(2R+2)
#
# cuya algebra de orbita es C_R, y de ahi que Phi_{2,R} sea un caracter SIMPLECTICO virtual: eso esta
# medido y es la base de toda la maquinaria.  Pero para t IMPAR el bloque es
#
#     mu_t = {1} u {xi^{+-j}}         ->   alfabeto = elemento de SO(2R'+1), SIN coset
#
# que es un elemento del grupo, no de la componente no trivial.  Por Littlewood, s_lambda restringido
# a SO(2R'+1) se descompone en caracteres ORTOGONALES, no simplecticos.  O sea que el borrador del
# paper afirma "caracter simplectico virtual" para las DOS paridades y puede estar mal en una.
#
# Esto NO se decide leyendo: se decide pelando en las dos bases y viendo cual cierra.
#
# LO QUE SE MIDE, sobre Phi_{1,R'} (el objeto impar sin especializar):
#   B_C   pelar en la base {sp_mu} de C_{R'}   -> ¿coeficientes ENTEROS?  ¿resto CERO?
#   B_B   pelar en la base {so_mu} de B_{R'}   -> idem
#
# CONTROLES
#   C0  ACEPTACION.  El mismo test sobre el objeto PAR Phi_{2,R}, donde ya sabemos la respuesta:
#       tiene que cerrar en C y NO tener por que cerrar en B.  Si el control no reprodujera lo
#       conocido, el instrumento no vale y no se dice nada del impar.
#   C1  no vacuidad: n impreso, y las dos columnas por separado.  Que una base cierre no dice nada
#       de la otra hasta que se prueban las dos.
#   C2  se prueba en VARIAS beta, no en una: una sola forma puede cerrar en las dos por casualidad.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage odd_basis.sage

import itertools, sys
from collections import defaultdict


def phi_de(beta, tt, nvar):
    N = tt + 2 * nvar
    K = CyclotomicField(tt) if tt > 2 else QQ
    zeta = K.gen() if tt > 2 else K(1 if tt == 1 else -1)
    L = LaurentPolynomialRing(K, nvar, 'z')
    zs = L.gens()
    x = [L(K(zeta) ** k) for k in range(tt)] + [g ** e for g in zs for e in (1, -1)]
    delta = list(range(N - 1, -1, -1))
    def alt(ex):
        return matrix(L, N, N, lambda i, j: x[i] ** ex[j]).determinant()
    den = alt(delta)
    if den == 0:
        return None
    q = alt(list(beta)) / den
    try:
        q = L(q)
    except Exception:
        return "NO-POL"
    return {tuple(e) if hasattr(e, '__iter__') else (e,): c
            for e, c in zip(q.exponents(), q.coefficients()) if c != 0}


_CH = {}
def car(mu, tipo, rk):
    """caracter irreducible de tipo 'C' o 'B' y rango rk, como dict sobre Z^rk."""
    key = (tuple(mu), tipo, rk)
    if key not in _CH:
        W = WeylCharacterRing("%s%d" % (tipo, rk))
        el = W(W.space().from_vector(vector(list(mu))))
        d = {}
        for wt, mult in el.weight_multiplicities().items():
            k = tuple(int(v) for v in wt.to_vector())
            d[k] = d.get(k, 0) + mult
        _CH[key] = d
    return _CH[key]


def pelar(P, tipo, rk, tope=3000):
    """P = sum c_mu chi_mu en la base dada.  Devuelve (coefs, resto)."""
    P = {e: QQ(c) for e, c in P.items() if c != 0}
    out = {}
    for _ in range(tope):
        if not P:
            return out, {}
        dom = [e for e in P if list(e) == sorted(e, reverse=True) and min(e) >= 0]
        if not dom:
            return out, P
        mu = max(dom, key=lambda e: (sum(e), e))
        c = P[mu]
        out[mu] = out.get(mu, 0) + c
        for k, v in car(mu, tipo, rk).items():
            nv = P.get(k, 0) - c * v
            if nv == 0:
                P.pop(k, None)
            else:
                P[k] = nv
    return out, P


def informe(P, rk, etiqueta):
    fila = []
    for tipo in ("C", "B"):
        co, resto = pelar(P, tipo, rk)
        co = {k: v for k, v in co.items() if v != 0}
        ent = all(QQ(v).denominator() == 1 for v in co.values())
        fila.append((len(co), len(resto), ent))
    print("   %-34s | C: %3d coefs, resto %3d, enteros %-3s | B: %3d coefs, resto %3d, enteros %s"
          % (etiqueta, fila[0][0], fila[0][1], "si" if fila[0][2] else "NO",
             fila[1][0], fila[1][1], "si" if fila[1][2] else "NO"))
    sys.stdout.flush()
    return fila


print("=" * 124)
print("¿EN QUE BASE VIVE EL OBJETO?   pelado en C_{rk} (simplectico) y en B_{rk} (ortogonal impar)")
print("=" * 124)
print("")
print("  C0  ACEPTACION -- el objeto PAR, donde la respuesta se conoce: tiene que cerrar en C")
print("")
print("   %-34s | %-38s | %s" % ("beta   (t=2, R=3, N=8)", "base C_3", "base B_3"))
print("   " + "-" * 118)
R = 3
for b in [(10, 9, 7, 4, 3, 2, 1, 0), (12, 10, 8, 5, 3, 2, 1, 0), (13, 11, 8, 6, 4, 2, 1, 0)]:
    P = phi_de(b, 2, R)
    if P in (None, "NO-POL") or not P:
        print("   %-34s | ---" % str(b))
        continue
    informe(P, R, str(b))

print("")
print("  N1  EL OBJETO IMPAR,  Phi_{1,R'}  con R' = 3,  N = 1 + 2R' = 7")
print("")
print("   %-34s | %-38s | %s" % ("beta   (t=1, R'=3, N=7)", "base C_3", "base B_3"))
print("   " + "-" * 118)
Rp = 3
for b in [(9, 7, 5, 3, 2, 1, 0), (11, 8, 6, 4, 3, 1, 0), (10, 9, 6, 4, 2, 1, 0),
          (12, 10, 7, 5, 3, 2, 0), (8, 6, 5, 3, 2, 1, 0)]:
    P = phi_de(b, 1, Rp)
    if P in (None, "NO-POL") or not P:
        print("   %-34s | ---  (nulo o no polinomio)" % str(b))
        continue
    informe(P, Rp, str(b))

print("")
print("=" * 124)
print("  LECTURA, escrita ANTES de correr:")
print("    * si el par cierra en C y el impar cierra en B  ->  el borrador del paper esta MAL en el")
print("      impar y hay que reescribir la seccion: la mitad impar ramifica por el lado ORTOGONAL.")
print("    * si los dos cierran en C  ->  el borrador vale y el impar es de verdad el companion.")
print("    * si el impar cierra en las DOS, no dice nada: hay que mirar cual da coeficientes enteros")
print("      y resto cero de forma estable sobre varias beta, que es la columna que se imprime.")
print("=" * 124)
print("DONE")
