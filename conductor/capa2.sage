# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# capa2.sage -- the second layer: is it the square of the first?
# What: (A) p not dividing N: W_1, W_2, and dim(V*V) with V = gr_1(A) as a subspace of F_p[y]/(Phi_{q'});
#       (B) the ladder zone p | N for p odd: profiles above P = p^{v_p(N)}.
# Run: sh sg.sh capa2.sage
load("inst.sage")

def capas_sub(I, K):
    """gr_i(A) for i < K, each as a matrix over F_p whose rows span it inside F_p[y]/(Phi)."""
    R, F, phin = I.R, I.F, I.phin
    rows = [R(list(r)) for r in I.B.rows()]
    dphi = phin.degree()
    out = []
    for i in range(K):
        Pi = phin**i
        if i == 0:
            elems = rows
        else:
            M = matrix(F, [I.vec(f % Pi, Pi.degree()) for f in rows])
            ker = M.left_kernel().basis_matrix()
            elems = [sum(ker[k][j] * rows[j] for j in range(len(rows))) for k in range(ker.nrows())]
        if not elems:
            out.append(matrix(F, 0, dphi)); continue
        img = matrix(F, [I.vec((f // Pi) % phin, dphi) for f in elems]).echelon_form()
        out.append(img.matrix_from_rows(range(img.rank())))
    return out

def producto(I, V):
    """span of the pairwise products of a basis of V, inside F_p[y]/(Phi)."""
    R, F, phin = I.R, I.F, I.phin
    dphi = phin.degree()
    b = [R(list(r)) for r in V.rows()]
    pr = [I.vec((b[i] * b[j]) % phin, dphi) for i in range(len(b)) for j in range(i, len(b))]
    return matrix(F, pr).rank() if pr else 0

def hs(n, c):
    c %= n
    if 2 * c == n: return 0
    return c if 2 * c < n else c - n

def familia_mult(m, n):
    """(N o B, centrada) por cardinal x punto medio"""
    S = [0] * n; cnt = [0] * n
    for j in range(-m, m + 1):
        S[j % n] += j; cnt[j % n] += 1
    centrada = all(S[c] == cnt[c] * hs(n, c) for c in range(n) if cnt[c])
    return (cnt[1] if centrada else cnt[1] // 2), centrada

print("=== (A) p not dividing N: is the second layer the square of the first? ===")
print("m     q      p  q'  d   W_1  W_2  dim(V*V) | W_2=d  W_2=dim(V*V)")
tot = bad_full = bad_prod = 0
casos = []
for qp in [5, 7, 9, 11, 13, 15, 16, 17, 19, 20, 21, 23, 24, 25]:
    d = euler_phi(qp) // 2
    for p in [3, 5, 7]:
        if qp % p == 0: continue
        v = 1
        while euler_phi(p**v) < 3: v += 1
        ms = []
        for N in [1, 2, 3]:
            if N % p == 0: continue
            ms += [N * qp, N * qp - 1]
            if (N * qp) % 2: ms.append((N * qp - 1) // 2)
            if qp % 2 == 0: ms += [N * (qp // 2), N * (qp // 2) - 1]
        for m in sorted(set(x for x in ms if x >= 2)):
            Nm, cen = familia_mult(m, qp)
            if Nm % p == 0: continue
            h = 2 * m
            if h % qp and (h + 1) % qp and (h + 2) % qp: continue
            vv = v
            while 2 * m >= p**vv * qp: vv += 1
            casos.append((m, p**vv * qp, p, qp, d))
for (m, q, p, qp, d) in casos:
    I = Img(q, p, m, 3)
    W = I.W()
    G = capas_sub(I, 3)
    dv = producto(I, G[1])
    full = (W[2] == d); prod = (W[2] == dv)
    tot += 1; bad_full += (not full); bad_prod += (not prod)
    print("%-5d %-6d %-2d %-3d %-3d %-4d %-4d %-8d | %-6s %s" % (
        m, q, p, qp, d, W[1], W[2], dv, full, prod))
    sys.stdout.flush()
print("A TOTAL=%d ; W_2 < d en %d ; W_2 != dim(V*V) en %d" % (tot, bad_full, bad_prod))

print("")
print("=== (B) the ladder zone, p odd: the layers above P ===")
print("m     q     p  q'  d   P  E  | W (i < E)          | W_1(A_0) | above P")
for p in [3, 5, 7]:
    for qp in [5, 7, 11, 13, 17, 19, 23, 25, 27]:
        if qp % p == 0: continue
        for M in [x for x in [1, 2, 3, 4] if x % p]:
            N = p * M
            if (N * qp) % 2 == 0: continue
            m = (N * qp - 1) // 2
            P = p
            v = 2
            while 2 * m >= p**v * qp or euler_phi(p**v) < P + 2: v += 1
            q = p**v * qp
            E = euler_phi(p**v)
            K = min(E, P + 3)
            t0 = walltime()
            I = Img(q, p, m, K)
            W = I.W()
            m0 = (M * qp - 1) // 2; q0 = q // P
            W10, _ = moment_rank(qp, m0, p)
            d = euler_phi(qp) // 2
            # candidatos para W_{P+1}: rangos mod p de matrices de segundo momento
            def rank2(mm, con_constantes):
                S2 = [0] * qp
                for j in range(-mm, mm + 1):
                    S2[j % qp] += j**2
                un = [u for u in range(1, qp) if gcd(u, qp) == 1 and u < qp - u]
                filas = [[S2[(inverse_mod(u, qp) * c) % qp] for c in range(qp)] for u in un]
                if con_constantes:
                    filas.append([1] * qp)
                return matrix(GF(p), filas).rank()
            def rank_hat2(con_constantes):
                """modulo generado por la funcion par c -> hat(c)^2, sin la constante del segundo momento"""
                f = [hs(qp, c)**2 for c in range(qp)]
                un = [u for u in range(1, qp) if gcd(u, qp) == 1 and u < qp - u]
                filas = [[f[(inverse_mod(u, qp) * c) % qp] for c in range(qp)] for u in un]
                if con_constantes:
                    filas.append([1] * qp)
                return matrix(GF(p), filas).rank()
            r2r = rank2(m0, False)
            h2, h2c = rank_hat2(False), rank_hat2(True)
            WP1 = W[P + 1] if len(W) > P + 1 else None
            ok = "hat^2" if WP1 == h2 else ("hat^2+cte" if WP1 == h2c else ("S_2^0" if WP1 == r2r else "FALLA"))
            print("%-5d %-6d %-2d %-3d %-3d %-2d %-2d | %-24s | %-4d | W_{P+1}=%-3s  <S_2^0>=%-3d <hat^2>=%-3d +cte=%-3d  %s" % (
                m, q, p, qp, d, P, E, W, W10, WP1, r2r, h2, h2c, ok))
            sys.stdout.flush()
