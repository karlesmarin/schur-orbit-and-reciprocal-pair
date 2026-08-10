# Authors: Carles Marin + Claude (AI assistant).
# I was about to write into the paper that Karmakar's case (B) "vanishes at r=3".
# That was an artefact of MY range (|lambda| <= 9).  Case (B) at k=1 needs N-1 of
# the N beta numbers to share a parity, which forces |lambda| to grow with N.
# Find the SMALLEST lambda in case (B) for each r, exactly.
def caseB_min(N):
    # beta: N distinct nonnegative integers, N-1 of one parity, 1 of the other,
    # minimising sum(beta).  Take the N-1 smallest of one parity, plus the
    # smallest available of the other.
    best = None
    for par in (0, 1):
        same = [2*i + par for i in range(N-1)]
        other = [x for x in range(0, 4*N) if x % 2 != par and x not in same]
        b = sorted(same + [other[0]], reverse=True)
        rho = list(range(N-1, -1, -1))
        lam = [b[i] - rho[i] for i in range(N)]
        if all(lam[i] >= lam[i+1] for i in range(N-1)) and lam[-1] >= 0:
            cand = (sum(lam), tuple(x for x in lam if x), tuple(b))
            if best is None or cand[0] < best[0]: best = cand
    return best

print("  smallest lambda in Karmakar's case (B) at k = 1")
print("  %-5s %-5s %-10s %-28s %s" % ("r", "N", "min|lambda|", "lambda", "beta"))
for r in (1, 2, 3, 4, 5):
    N = 2*r+2
    m = caseB_min(N)
    print("  %-5d %-5d %-10d %-28s %s" % (r, N, m[0], str(list(m[1])), list(m[2])))
print("")
print("  So case (B) does NOT die -- it moves out of reach.  The honest statement")
print("  is that its smallest member grows, not that the case empties.")
