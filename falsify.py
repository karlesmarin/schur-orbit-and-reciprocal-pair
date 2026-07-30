"""Falsification controls for the general-m law. Authors: Carles Marin, Claude (AI assistant)."""
from mpmath import mp, mpf, mpc, exp, sinh, pi
from law_control import mydet, partitions, law_prediction
mp.dps = 40

def phi_alt(lam, m, theta, mode):
    N = m + 2
    t = exp(theta); w = exp(2j*pi/m)
    if mode == "orbit":   xs = [w**k for k in range(m)] + [t, 1/t]
    elif mode == "coset": z = exp(1j*pi/m); xs = [z*w**k for k in range(m)] + [t, 1/t]   # odd powers of a 2m-th root (Ayyer-Kumari 2501.00275 s3)
    elif mode == "free":  xs = [w**k for k in range(m)] + [t, exp(mpf("0.77"))]          # free pair, not reciprocal
    lam = list(lam) + [0]*(N-len(lam))
    beta = [lam[j] + N-1-j for j in range(N)]
    num = [[xs[i]**beta[j] for j in range(N)] for i in range(N)]
    den = [[xs[i]**(N-1-j) for j in range(N)] for i in range(N)]
    return mydet(num, N)/mydet(den, N)

theta = mpf("0.41")
for mode in ("orbit", "coset", "free"):
    ok = bad = 0
    for m in range(2, 7):
        for n in range(0, 11):
            for lam in partitions(n, m+2):
                phi = phi_alt(lam, m, theta, mode)
                pred = law_prediction(lam, m, theta)
                if pred is None or abs(pred[1]) < 1e-18:
                    good = abs(phi) < 1e-12
                else:
                    r = phi/pred[1]
                    good = abs(r.imag) < 1e-12 and abs(abs(r.real)-1) < 1e-12
                ok += good; bad += (not good)
    print(f"{mode:6s}: law holds {ok:5d}   law FAILS {bad:5d}")
