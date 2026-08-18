# -*- coding: utf-8 -*-
# Authors: Carles Marin, Claude (AI assistant).
# El que encontro el defecto del 2026-08-15: el enunciado impreso de lem:step contaba M solo sobre
# las clases de exceso.  Ver [[the-script-may-be-right-and-the-paper-wrong]].
"""Does lem:step hold with M read as the paper states it (x in g, excess classes only),
or only with M read as the script does (x in P, ALL frozen values, one per class)?

w(S) = (-1)^{sum_{j in S} j + inv(b_S)}, S the frozen columns.  We move the choice in ONE excess
class from u down to v and compare w(g')/w(g) against both readings.
"""
import itertools

def inv(word):
    return sum(1 for a in range(len(word)) for b in range(a + 1, len(word)) if word[a] > word[b])

def audit(t, r, W):
    N = t + 2 * r
    diff_cases = 0
    paper_wrong = script_wrong = 0
    witness = None
    for mid in itertools.combinations(range(1, W + 1), N - 1):
        beta = tuple(sorted(mid, reverse=True)) + (0,)
        cls = {}
        for idx, v in enumerate(beta):
            cls.setdefault(v % t, []).append((v, idx))
        if len(cls) < t:
            continue                      # occupancy hypothesis
        excess = [i for i in cls if len(cls[i]) >= 2]
        if not excess:
            continue
        # a choice g: one value per EXCESS class; frozen set P = g + all singleton values
        singles = [cls[i][0] for i in cls if len(cls[i]) == 1]
        for choice in itertools.product(*[cls[i] for i in excess]):
            g = list(choice)
            P = g + singles
            Scols = sorted(idx for (_, idx) in P)
            bS = [beta[j] % t for j in Scols]
            w = (-1) ** (sum(j + 1 for j in Scols) + inv(bS))
            for ci, i in enumerate(excess):
                for (u, uidx) in cls[i]:
                    for (v, vidx) in cls[i]:
                        if not (u > v) or (u, uidx) != choice[ci]:
                            continue
                        g2 = list(choice)
                        g2[ci] = (v, vidx)
                        P2 = list(g2) + singles
                        Scols2 = sorted(idx for (_, idx) in P2)
                        bS2 = [beta[j] % t for j in Scols2]
                        w2 = (-1) ** (sum(j + 1 for j in Scols2) + inv(bS2))
                        truth = w2 * w        # ratio, both +-1

                        B = sum(1 for x in beta if v < x < u)
                        M_paper = sum(1 for (x, _) in g if v < x < u)
                        M_script = sum(1 for (x, _) in P if v < x < u)
                        pred_p = (-1) ** (1 + B + M_paper)
                        pred_s = (-1) ** (1 + B + M_script)
                        if M_paper != M_script:
                            diff_cases += 1
                            if pred_p != truth:
                                paper_wrong += 1
                                if witness is None:
                                    witness = (beta, i, u, v, B, M_paper, M_script, truth, pred_p)
                            if pred_s != truth:
                                script_wrong += 1
    return diff_cases, paper_wrong, script_wrong, witness


print("  t  r  W | casos donde las dos lecturas DIFIEREN | falla 'x in g' | falla 'x in P'")
for (t, r, W) in [(3, 1, 9), (4, 1, 10), (5, 1, 11), (3, 2, 11), (4, 2, 12), (5, 2, 12), (6, 2, 13)]:
    d, pw, sw, wit = audit(t, r, W)
    print("  %2d %2d %2d | %8d | %8d | %8d" % (t, r, W, d, pw, sw))
    if wit:
        print("       testigo: beta=%s clase=%d u=%d v=%d  B=%d  M(g)=%d M(P)=%d  real=%+d  'x in g' da %+d"
              % (wit[0], wit[1], wit[2], wit[3], wit[4], wit[5], wit[6], wit[7], wit[8]))
