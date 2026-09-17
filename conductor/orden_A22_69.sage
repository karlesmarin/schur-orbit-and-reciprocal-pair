# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# What: A_22(69) = Z[e_1..e_22(eta_1..eta_22)] inside O = Z[zeta_69+zeta_69^-1]: index, conductor, layers W_i at 3, CM type (two ways).
# Run (in this directory): sage orden_A22_69.sage
import time
load("orden_generado.sage")
t0 = time.time()
q = 69; m = 22
C = CyclotomicField(q); z = C.gen()
fpol = (z + 1/z).minpoly()
L = NumberField(fpol, 'a'); a = L.gen()
d = L.degree()
print("deg K =", d, " minpoly disc factor:", factor(fpol.discriminant()))
OK = L.maximal_order()
Za = L.order(a)
print("[O_K : Z[a]] =", Za.index_in(OK))

eta = [L(2), a]
for j in range(1, m):
    eta.append(a*eta[j] - eta[j-1])
eta = eta[1:m+1]  # eta_1..eta_m
R = PolynomialRing(L, 'X'); X = R.gen()
P = prod(X - e for e in eta)
e = [(-1)**k * P[m-k] for k in range(1, m+1)]
assert all(ei.is_integral() for ei in e)

def co(x): return vector(QQ, list(L(x)))
def el(v): return L(list(v))
def hnfQ(vs):
    M = matrix(QQ, [list(v) for v in vs])
    N = lcm([t.denominator() for t in M.list()])
    H = matrix(ZZ, N*M).hermite_form(include_zero_rows=False)
    return H / N
def idx(H):  # index in O = Z^d of a full-rank lattice with basis rows H
    assert H.nrows() == d
    return abs(H.det())

# --- A by lattice closure (own loop) and by og_orden
H = hnfQ([co(1)])
while True:
    base = [el(r) for r in H.rows()]
    H2 = hnfQ([co(b) for b in base] + [co(g*b) for g in e for b in base])
    if H2 == H: break
    H = H2
HA = H
iA = idx(HA)
print("[O:A] (closure) =", iA, "=", factor(iA))
Aord = og_orden(L, e)
print("[O_K:A] (og_orden, index_in maximal order) =", Aord.index_in(OK))

def kernel_mod(Cm, D):
    # {x in Z^n : x*Cm == 0 mod D}, Cm integer n x c
    n, c = Cm.nrows(), Cm.ncols()
    big = block_matrix(ZZ, [[Cm, identity_matrix(ZZ, n)], [D*identity_matrix(ZZ, c), zero_matrix(ZZ, c, n)]])
    Hb = big.hermite_form()
    rows = [r[c:] for r in Hb.rows() if r[:c] == 0 and r[c:] != 0]
    return matrix(ZZ, rows).hermite_form(include_zero_rows=False)

def sub_lattice(Hbig_rows_list, Htarget, Hstart=None):
    # {x in lattice(Hstart) (default Z^d) : x*Bk in lattice(Htarget) for all k}, Bk rational d x d.
    # Sequential intersection: F <- {x in F : x*Bk in target}, each a (d+d)x(d+d) HNF.
    Hinv = Htarget.inverse()
    F = identity_matrix(QQ, d) if Hstart is None else matrix(QQ, Hstart)
    for B in Hbig_rows_list:
        Cm = F*B*Hinv
        D = lcm([t.denominator() for t in Cm.list()])
        Y = kernel_mod(matrix(ZZ, D*Cm), D)
        F = hnfQ(list(matrix(QQ, Y)*F))
    return F

def multmat(x):  # row-convention matrix of y -> y*x on O coords
    return matrix(QQ, [co(a**i * x) for i in range(d)])
def intersect(H1, H2):
    Y = sub_lattice([H1], H2)  # y with y*H1 in lattice(H2)
    return (Y*H1).hermite_form() if H1.base_ring() == ZZ else hnfQ(list(Y*H1))
def lat_sum(*Hs):
    return hnfQ([r for Hh in Hs for r in Hh.rows()])
def lat_prod(H1, H2):
    return hnfQ([co(el(r)*el(s)) for r in H1.rows() for s in H2.rows()])
I = identity_matrix(QQ, d)

# --- conductor f = {x in O : x a^k in A, 0<=k<d}
Ta = multmat(a)
Hf = sub_lattice([Ta**k for k in range(d)], HA)
Hf = matrix(QQ, Hf)
print("[O:f] =", idx(Hf), "=", factor(idx(Hf)))
print("f == 3O :", Hf == hnfQ(list(3*I)))
print("3O subset A :", all(x in ZZ for x in (3*I*HA.inverse()).list()))

# --- radical r above 3
fb = fpol.change_ring(GF(3)).factor()
print("minpoly mod 3 =", [(g.degree(), ex) for g, ex in fb])
assert len(fb) == 1 and fb[0][1] == 2
g = fb[0][0].change_ring(ZZ)
ga = g(a)
Hr = hnfQ(list(3*I) + [co(a**i * ga) for i in range(d)])
print("[O:r] =", factor(idx(Hr)))
Hr2 = lat_prod(Hr, Hr)
print("r^2 == 3O :", Hr2 == hnfQ(list(3*I)))
rpow = [hnfQ(list(I)), Hr, Hr2, lat_prod(Hr2, Hr), lat_prod(Hr2, Hr2)]
e_exp = min(i for i in range(len(rpow)) if all(x in ZZ for x in (rpow[i]*HA.inverse()).list()))
print("smallest e with r^e subset A:", e_exp, "; f == r^e :", Hf == rpow[e_exp])
W = []
for i in range(4):
    Ai = intersect(HA, rpow[i])
    S = lat_sum(Ai, rpow[i+1])
    w = ZZ(idx(rpow[i+1]) / idx(S)).log(3)
    assert 3**w == idx(rpow[i+1]) / idx(S)
    W.append(w)
print("layers W_0..W_3 =", W)

# --- CM type via colon (A :_O m), m = A cap r
Hm = intersect(HA, Hr)
print("[A:m] =", factor(idx(Hm)/iA))
Hcol = sub_lattice([multmat(el(r)) for r in Hm.rows()], HA)
Hcol = matrix(QQ, Hcol)
print("A subset (A:_O m):", all(x in ZZ for x in (HA*Hcol.inverse()).list()))
colidx = iA / idx(Hcol)
print("[(A:_O m):A] =", colidx, "; log_3 =", ZZ(colidx).log(3), "; |A/m| = 3^%d" % W[0])

# --- independent: canonical module omega = trace dual of A; type = mu(omega) = log_|k| [omega : m omega]
Bel = [el(r) for r in HA.rows()]
G = matrix(QQ, d, d, lambda i, j: (Bel[i]*Bel[j]).trace())
Hom = hnfQ(list(G.inverse()*HA))
Hmom = lat_prod(Hm, Hom)
ratio = idx(Hmom)/idx(Hom)
print("[omega : m omega] =", ratio, "; log_3 =", ZZ(ratio).log(3))
print("Gorenstein (colon) :", ZZ(colidx).log(3) == W[0], " Gorenstein (mu(omega)=1):", ZZ(ratio).log(3) == W[0])
print("time %.1fs" % (time.time()-t0))
