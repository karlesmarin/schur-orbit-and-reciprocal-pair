# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# cruce.sage -- crossing the formulas of the conductor note, reduced to minimal form.
# What: (A) one identity behind both mechanisms: 2*hat_s_n = s_n o [2] - s_n for EVERY n, and for n even
#           it equals the translation s_n(c+n/2); (B) the first moment is count x midpoint on the four
#           stable families; (C) the two branches (T_l - 1) and (T_l - l) of the same operator family:
#           ranks over Q and mod p against the character count.   Run: sh sg.sh cruce.sage

hm = {}
for ln in open("stickelberger_smith_OUT.txt"):
    t = ln.split()
    if t and t[0].isdigit():
        hm[int(t[0])] = Integer(t[4])

def s(n, c):
    c %= n
    return 0 if c == 0 else 2 * c - n

def hs(n, c):
    c %= n
    if 2 * c == n: return 0
    return c if 2 * c < n else c - n

def reps(n):
    return [u for u in range(1, n) if gcd(u, n) == 1 and u < n - u]

def mat(n, f, R=ZZ):
    return matrix(R, [[f(inverse_mod(u, n) * c) for c in range(n)] for u in reps(n)])

print("=== (A) 2*hat_s_n(c) = s_n(2c) - s_n(c) for every n; = s_n(c+n/2) for n even ===")
bad = badpar = 0
for n in range(3, 101):
    for c in range(n):
        if 2 * hs(n, c) != s(n, 2 * c) - s(n, c):
            bad += 1
        if n % 2 == 0 and 2 * hs(n, c) != s(n, c + n // 2):
            badpar += 1
print("n=3..100: fallos de la identidad general = %d ; fallos de la version par = %d" % (bad, badpar))
print("corolario: para n par, s_n(2c) - s_n(c) = s_n(c+n/2) (el operador [2] degenera en una traslacion)")
print("")

print("=== (B) first moment = count x midpoint, on the four stable families ===")
print("family            q'  m     | N_c const | midpoint = | S = N_c*mu")
tot = bad = 0
for qp in [7, 8, 9, 12, 13, 15, 16, 20, 21, 24, 25]:
    for (fam, ms) in [("q'|h      m=Nq'", [qp, 2 * qp, 3 * qp]),
                      ("q'|h+2    m=Nq'-1", [qp - 1, 2 * qp - 1, 3 * qp - 1]),
                      ("q'|h+1    2m+1=Nq'", [(qp - 1) // 2, (3 * qp - 1) // 2] if qp % 2 else []),
                      ("half per. m=Br", [qp // 2, 3 * (qp // 2)] if qp % 2 == 0 else []),
                      ("half per. m=Br-1", [qp // 2 - 1, 3 * (qp // 2) - 1] if qp % 2 == 0 else [])]:
        for m in ms:
            if m < 1: continue
            S = [0] * qp; N = [0] * qp
            for j in range(-m, m + 1):
                S[j % qp] += j; N[j % qp] += 1
            # midpoint of each class-progression
            mu = [(S[c] / N[c] if N[c] else None) for c in range(qp)]
            okprod = all(N[c] * mu[c] == S[c] for c in range(qp) if N[c])
            # which model: centred (mu = hat c) or off-centre (mu = s/2)
            centred = all(mu[c] == hs(qp, c) for c in range(qp) if N[c])
            offc = all(2 * mu[c] == s(qp, c) for c in range(qp) if N[c])
            Nset = sorted(set(N[c] for c in range(1, qp) if 2 * c != qp))
            tot += 1; bad += not (okprod and (centred or offc) and len(Nset) == 1)
            print("%-18s %2d %4d | %-9s | %-10s | %s" % (
                fam, qp, m, str(Nset), "hat c" if centred else ("s/2" if offc else "OTRA"),
                "OK" if okprod else "FALLA"))
print("FAMILIES TOTAL=%d fail=%d" % (tot, bad))
print("")

print("=== (C) the two branches of the same operator family, on odd q' ===")
print("S_l = (s(lc)-s(c))/2  [T_l - 1]   vs   C_l = (s(lc)-l*s(c))/(2q')  [T_l - l]")
print("q'  l  d  t=ord(l) -1 in <l> | rank_Q S_l  pred d-delta | rank_Q C_l  vals(C_l) | p-sweep")
tot = bad = 0
for qp in [5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35]:
    d = euler_phi(qp) // 2
    for l in [2, 3, 5, 7]:
        if gcd(l, qp) > 1: continue
        t = Mod(l, qp).multiplicative_order()
        inv = (Mod(-1, qp) in [Mod(l, qp)**k for k in range(t)])
        delta = 0 if inv else d / t
        # (T_l - l) is integral only for odd l; scale by 2 so that both branches are integer matrices
        sc = 1 if l % 2 else 2
        Sl = mat(qp, lambda c: (s(qp, l * c) - s(qp, c)) / 2)
        Cl = mat(qp, lambda c: sc * (s(qp, l * c) - l * s(qp, c)) / (2 * qp))
        rS = Sl.rank(); rC = Cl.rank()
        vals = sorted(set(Cl.list()))
        # mod p sweep, p odd, p not dividing q' * h^-
        msg = []
        for p in [3, 5, 7, 11, 13]:
            if qp % p == 0 or hm[qp] % p == 0: continue
            rSp = Sl.change_ring(GF(p)).rank(); rCp = Cl.change_ring(GF(p)).rank()
            # character prediction: number of odd chi with (chibar(l)-1) or (chibar(l)-l) nonzero mod p
            # orbit form: prod over an orbit of length t is (l^t - 1) for the second branch
            dropC = 0 if (l**t - 1) % p else "?"
            msg.append("p=%d:S=%d,C=%d%s" % (p, rSp, rCp, "" if rCp == d else ("(baja;p|l^t-1:%s)" % ((l**t - 1) % p == 0))))
        # (T_l - l)B_1 = (l-1)/2 - floor(l x): a step function with l levels, symmetric around 0
        niv = set(range(-((l - 1) // 2), (l - 1) // 2 + 1)) if l % 2 else set([-1, 0, 1])
        ok = (rS == d - delta) and (rC == d) and set(vals) <= niv
        tot += 1; bad += not ok
        print("%2d %2d %2d   %2d      %-5s | %2d  %2d | %2d  %-12s | %s %s" % (
            qp, l, d, t, inv, rS, d - delta, rC, str(vals), " ".join(msg), "OK" if ok else "FALLA"))
print("OPERATOR TOTAL=%d fail=%d" % (tot, bad))
print("")

print("=== (C2) control: the rank d - delta_l needs p ODD ===")
print("the formula is stated over Q and modulo odd p not dividing n*h^-; at p=2 it fails")
mal2 = 0
for (n, l) in [(5, 4), (5, 2), (7, 2), (9, 2), (11, 2), (13, 4), (7, 4)]:
    if gcd(l, n) > 1: continue
    d = euler_phi(n) // 2
    t = Mod(l, n).multiplicative_order()
    inv = (Mod(-1, n) in [Mod(l, n)**k for k in range(t)])
    delta = 0 if inv else d / t
    Sl = mat(n, lambda c: (s(n, l * c) - s(n, c)) / 2)
    r2 = Sl.change_ring(GF(2)).rank()
    coincide = (r2 == d - delta)
    mal2 += (not coincide)
    print("n=%2d l=%d d=%2d | prediccion d-delta=%2s | rango mod 2 = %d  %s%s" % (
        n, l, d, d - delta, r2, "coincide" if coincide else "FALLA (por eso p impar)",
        "" if n % 2 == 0 else "   [h^-(%d)=%s]" % (n, hm.get(n, "?"))))
print("controles en p=2: %d de %d fallan la formula, con p no dividiendo n*h^-" % (mal2, 7))
print("")

print("=== (D) second order: which moment does the next layer see? ===")
print("eta_j = gamma_j + j*delta_j*pi + [ (j^2/2) gamma_j - (j/2) delta_j ] pi^2 + O(pi^3)")
tot = bad = 0
for n in [5, 7, 9, 11, 13]:
    K = CyclotomicField(n); w = K.gen()
    R.<pi> = PowerSeriesRing(K, default_prec=4)
    rho = 1 + pi
    for j in range(1, 3 * n):
        eta = w**j * rho**j + w**(-j) * rho**(-j)
        g = w**j + w**(-j); dl = w**j - w**(-j)
        pred = g + j * dl * pi + (j**2 / 2 * g - j / 2 * dl) * pi**2
        tot += 1; bad += ((eta - pred).valuation() < 3)
print("expansion checked on %d pairs (n,j); mismatches = %d" % (tot, bad))
print("consequence: the pi-coefficient carries the FIRST moment (odd, weight delta),")
print("the pi^2-coefficient carries the SECOND moment (even, weight gamma) and the first again.")
print("")
print("second moment on the centred families: S2(c) = sum_{|j|<=m, j=c} j^2")
for (n, m) in [(5, 7), (5, 12), (7, 10), (7, 17), (9, 13), (11, 16), (13, 19)]:
    if (2 * m + 1) % n: continue
    N = (2 * m + 1) // n
    S2 = [0] * n
    for j in range(-m, m + 1):
        S2[j % n] += j**2
    base = S2[0]
    ok = all(S2[c] == N * hs(n, c)**2 + base for c in range(n))
    print("q'=%2d m=%2d N=%d | S2(c) = N*hat_c^2 + S2(0) : %s  (S2(0)=%d)" % (n, m, N, ok, base))
