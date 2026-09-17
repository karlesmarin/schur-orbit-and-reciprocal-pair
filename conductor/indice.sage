# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# indice.sage -- closed formula for the index: v_p([O:A]) = 2d-1-W_1 whenever W_1 > d/2,
#                and its explicit form d-1+k_p (elementary divisors) or d-1+delta (factor at 2),
#                against the lattice indices of datos_indice.json.
# What: W_1 by moment_rank (inst.sage); k_p by the Smith form of M_{q'}; h^- from
#       stickelberger_smith_OUT.txt.   Run: sh sg.sh indice.sage
load("inst.sage")
import json

hm = {}
for ln in open("stickelberger_smith_OUT.txt"):
    t = ln.split()
    if t and t[0].isdigit():
        hm[int(t[0])] = Integer(t[4])

def s(n, c):
    c %= n
    return 0 if c == 0 else 2 * c - n

def reps(n):
    return [u for u in range(1, n) if gcd(u, n) == 1 and u < n - u]

def Mn(n):
    return matrix(ZZ, [[s(n, inverse_mod(u, n) * c) for c in range(n)] for u in reps(n)])

ED = {}
def ed_div(n, p):
    if n not in ED:
        ED[n] = [Integer(x) for x in Mn(n).elementary_divisors() if x != 0]
    return len([x for x in ED[n] if x % p == 0])

celdas = json.load(open("datos_indice.json"))["celdas"]

print("=== closed index formula against %d lattice cells ===" % len(celdas))
print("m   q    p  q'  d  W_1  | medido  2d-1-W_1 | familia     k_p/delta  explicito | ")
tot = bad = fuera = 0
fuera_list = []
expl_tot = expl_bad = 0
for cel in celdas:
    m, q = cel["m"], cel["q"]
    for pr in cel["primos"]:
        p, a, qp = pr["p"], pr["a"], pr["qprima"]
        if p == 2 or qp in (1, 2, 3, 4, 6):
            continue
        if pr["familia"] is None:      # outside the support: Theorem 3(c) gives a = 0
            if a != 0:
                print("SOPORTE FALLA m=%d q=%d p=%d a=%d" % (m, q, p, a)); bad += 1
            continue
        d = euler_phi(qp) // 2
        W1, _ = moment_rank(qp, m, p)
        # which of the two midpoints: centred (S = N*hat c) or off-centre (S = N*s/2)
        S = [0] * qp; cnt = [0] * qp
        for j in range(-m, m + 1):
            S[j % qp] += j; cnt[j % qp] += 1
        r = qp // 2 if qp % 2 == 0 else None
        centrada = all(S[c] == cnt[c] * (c if 2 * c < qp else (0 if 2 * c == qp else c - qp))
                       for c in range(qp) if cnt[c])
        if 2 * W1 <= d:
            fuera += 1; fuera_list.append((m, q, p, qp, d, W1, a)); continue
        pred = 2 * d - 1 - W1
        ok = (pred == a)
        tot += 1; bad += (not ok)
        # explicit form
        expl = "-"
        if centrada and qp % 2:                      # family q' | h+1: factor at 2
            N = (2 * m + 1) // qp
            if N % p and hm.get(qp, 0) % p:
                t = Mod(2, qp).multiplicative_order()
                inv = (Mod(-1, qp) in [Mod(2, qp)**k for k in range(t)])
                delta = 0 if inv else d / t
                expl = "delta=%s" % delta
                expl_tot += 1; expl_bad += (a != d - 1 + delta)
        else:                                        # sawtooth or half-period
            B = cnt[1] if not centrada else cnt[1]
            if B % p:
                k = ed_div(qp, p)
                expl = "k=%d" % k
                expl_tot += 1; expl_bad += (a != d - 1 + k)
        print("%-3d %-4d %-2d %-3d %-2d %-3d  | %-6d %-8d | %-11s %-10s %s" % (
            m, q, p, qp, d, W1, a, pred, "centrada" if centrada else "diente", expl,
            "OK" if ok else "FALLA"))
        sys.stdout.flush()
print("")
print("CERRADA: %d pares (celda,primo) con W_1 > d/2, fallos = %d" % (tot, bad))
print("EXPLICITA: %d con p no dividiendo N (o B) y, en el caso centrado, tampoco h^-; fallos = %d" % (expl_tot, expl_bad))
print("FUERA (W_1 <= d/2, la zona p|N y las caidas grandes): %d" % fuera)
for x in fuera_list:
    print("   fuera: m=%d q=%d p=%d q'=%d d=%d W_1=%d v_p=%d" % x)
