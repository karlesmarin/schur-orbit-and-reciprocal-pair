# -*- coding: utf-8 -*-
# LAS CUATRO NORMALIZACIONES, a la vez.  Emparejar con corraiz o con raiz, y tomar d = orden en el
# ADJUNTO o d = orden en G.  El G2 de NPP §7 no distingue el primer eje (d=2, factor 3 impar) y
# no tiene el segundo (alli el orden en G y en el adjunto coinciden).  Aqui se ven los cuatro.
import importlib.util
import itertools
import sys

sys.path.insert(0, r"E:\proyectos\Curiosity\research\orbit-pair\gates")
spec = importlib.util.spec_from_file_location(
    "q", r"E:\proyectos\Curiosity\research\orbit-pair\gates\npp_question81_v2.py")
# no ejecutamos el modulo entero (imprime); copiamos lo minimo
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

exec(open(r"E:\proyectos\Curiosity\research\orbit-pair\gates\npp_question81_v2.py",
          encoding="utf-8").read().split('print("=" * 100)')[0])

LMAX = 9
print("=" * 96)
print("LAS CUATRO NORMALIZACIONES DE LA QUESTION 8.1 SOBRE NUESTRA FAMILIA")
print("=" * 96)
print("   A = emparejar con la CORRAIZ de C_m (= raiz del dual B_m)")
print("   B = emparejar con la RAIZ de C_m (el propio grupo)")
print("   d_ad = orden en el grupo adjunto     d_G = orden en G")
print("")

for m in (1, 2, 3):
    t = 2 * m + 2
    phi = ciclotomico(t)
    a = list(range(1, m + 1))
    d_ad = None
    for k in range(1, t + 1):
        ex = [(k * ai) % t for ai in a]
        if all(e == 0 for e in ex) or (t % 2 == 0 and all(e == t // 2 for e in ex)):
            d_ad = k
            break
    d_G = t

    lams = [l for l in itertools.product(range(0, LMAX + 1), repeat=m)
            if all(l[i] >= l[i + 1] for i in range(m - 1))]
    res = {}
    n = 0
    for lam in lams:
        v = sp_valor(list(lam), m, a, t, phi)
        if v is None or isinstance(v, str):
            continue
        n += 1
        lr = [lam[i] + (m - i) for i in range(m)]
        for lect in ("A", "B"):
            for nomd, dd in (("d_ad", d_ad), ("d_G", d_G)):
                ok = (hits(lr, m, dd, lect) > 0) == (v == 0)
                res[(lect, nomd)] = res.get((lect, nomd), 0) + (1 if ok else 0)
    print("  m=%d  Sp(%d)  t=%d  d_ad=%d  d_G=%d   %d particiones"
          % (m, 2 * m, t, d_ad, d_G, n))
    for k in (("A", "d_ad"), ("A", "d_G"), ("B", "d_ad"), ("B", "d_G")):
        marca = "  <== EQUIVALENCIA" if res[k] == n else ""
        print("      %s con %-5s : %3d de %3d%s" % (k[0], k[1], res[k], n, marca))
    print("")
