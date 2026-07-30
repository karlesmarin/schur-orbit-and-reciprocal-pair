# -*- coding: utf-8 -*-
"""Two counts the paper's verification table quotes but no script printed.

TEST A -- Proposition on rectangles, read LITERALLY, over t=2, 1<=r<=4, 1<=c<=60 (240 cases):
          (i)   c even : one d equals 2; the other two equal c+2 for r<=3, all three equal 2 for r=4;
                         count = ((c+2)/2)^2  (resp. 1)
          (ii)  c odd, r=2 : count 0
          (iii) c odd, r!=2 : the other two are c+1 and c+3, count = +-(c+1)(c+3)/4
          The signed count is the z->1 value of the closed form, computed as a limit.

TEST B -- the four box families of the Example, c=1..4 (16 cases):
          1x3xc: floor((c+2)^2/4)   2x2xc: (c/2+1)^2 (c even), 0 (c odd)
          3x1xc: (-1)^c floor((c+2)^2/4)   4x0xc: (-1)^c

Authors: Carles Marin, Claude (AI assistant)."""
from mpmath import mp, mpf, nint
from theorem_full import closed_form, setup
mp.dps = 60
TH = mpf('1e-15')          # z -> 1; the closed form is a ratio of sinh's, no cancellation


def signed_count(r, c):
    """the z->1 value of the closed form for lambda = (c^r) at t=2."""
    lam = tuple([c] * r)
    v = closed_form(lam, 2, TH)
    if v is None:
        return None
    v = mp.mpmathify(v)
    if abs(mp.im(v)) > mpf('1e-20'):
        raise ValueError("non-real closed form at r=%d c=%d : %s" % (len(lam), lam[0], v))
    return int(nint(mp.re(v)))


def dtriple(r, c):
    """the multiset {d1,d2,d3} for lambda = (c^r) at t=2, or None if the shape is skipped."""
    from theorem_full import setup as _s
    st = _s(tuple([c] * r), 2)
    if st is None:
        return None
    beta, Ac, Bc = st
    a1, a2 = beta[Ac[0]], beta[Ac[1]]
    b1, b2 = beta[Bc[0]], beta[Bc[1]]
    return tuple(sorted([abs(a1 - a2), abs(b1 - b2), abs(a1 + a2 - b1 - b2)]))


# ---------------------------------------------------------------- TEST A
print("TEST A -- the rectangle proposition read literally, t=2, r<=4, c<=60")
okA = badA = skipA = 0
fails = []
for r in range(1, 5):
    for c in range(1, 61):
        got = signed_count(r, c)
        if got is None:
            skipA += 1
            continue
        if c % 2 == 0:
            pred = 1 if r == 4 else ((c + 2) // 2) ** 2
        elif r == 2:
            pred = 0
        else:                                    # case (iii), literally r != 2
            pred = (c + 1) * (c + 3) // 4
        if abs(got) == abs(pred):
            okA += 1
        else:
            badA += 1
            if len(fails) < 10:
                fails.append((r, c, got, pred, dtriple(r, c)))
print("   cases %d ; agree (up to sign) %d ; DISAGREE %d ; skipped %d"
      % (okA + badA, okA, badA, skipA))
for r, c, got, pred, d in fails:
    print("      r=%d c=%-3d closed form = %-6d  proposition predicts %-6d  d-triple=%s"
          % (r, c, got, pred, str(d)))

# ---------------------------------------------------------------- TEST B
print("\nTEST B -- the four box families of the Example, c=1..4")
okB = badB = 0
for c in range(1, 5):
    for r, name, pred in [
            (1, "1x3xc", ((c + 2) ** 2) // 4),
            (2, "2x2xc", ((c // 2 + 1) ** 2 if c % 2 == 0 else 0)),
            (3, "3x1xc", ((-1) ** c) * (((c + 2) ** 2) // 4)),
            (4, "4x0xc", (-1) ** c)]:
        got = signed_count(r, c)
        tag = "ok" if (got is not None and got == pred) else "MISMATCH"
        if tag == "ok":
            okB += 1
        else:
            badB += 1
            print("      %s c=%d : closed form = %s , example says %d" % (name, c, got, pred))
print("   boxes checked %d ; agree %d ; MISMATCH %d" % (okB + badB, okB, badB))
