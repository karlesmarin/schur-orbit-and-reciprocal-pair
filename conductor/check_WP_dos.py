"""Check of the corollary on W_P (label cor:WPdos): W_P = d if -1 lies in <2> mod q', and d - d/t otherwise, when p does not
divide h^-(Q(zeta_q')).

Reads descent_OUT_all.txt (the measured W_P of the descent) and stickelberger_smith_OUT.txt (exact h^-). They are
looked for next to this script (flat archive layout) and otherwise in ../../anc (the working tree). Cases with
p | h^- are counted apart; p | M does not occur, since the descent requires p not dividing M.
"""
import re
from pathlib import Path

AQUI = Path(__file__).resolve().parent


def fichero(nombre):
    for d in (AQUI, AQUI.parents[1] / "anc"):
        if (d / nombre).exists():
            return d / nombre
    raise SystemExit("not found: " + nombre)


hm = {}
for linea in fichero("stickelberger_smith_OUT.txt").read_text(encoding="utf-8").splitlines():
    t = linea.split()
    if len(t) > 4 and t[0].isdigit():
        hm[int(t[0])] = int(t[4])


def potencias2(n):
    x, vistas = 1, []
    while True:
        vistas.append(x)
        x = x * 2 % n
        if x == 1:
            return vistas


pat = re.compile(r"p=\s*(\d+) a=\s*(\d+) n=\s*(\d+) M=\s*(\d+) v=\s*(\d+) .*?d=\s*(\d+) \|.*?W_P=(\d+)")
vistos, ok, fallos, sin_h = set(), 0, [], 0
fuera_total, fuera_difiere = 0, 0
for linea in fichero("descent_OUT_all.txt").read_text(encoding="utf-8").splitlines():
    r = pat.search(linea)
    if not r:
        continue
    p, a, n, M, v, d, WP = map(int, r.groups())
    clave = (p, a, n, M, v)
    if clave in vistos:
        continue
    vistos.add(clave)
    if n not in hm:
        sin_h += 1
        continue
    pot = potencias2(n)
    pred = d if (n - 1) in pot else d - d // len(pot)
    if M % p == 0 or hm[n] % p == 0:
        fuera_total += 1
        fuera_difiere += pred != WP
        continue
    if pred == WP:
        ok += 1
    else:
        fallos.append(linea)
print(f"casos distintos {len(vistos)}; en hipotesis {ok + len(fallos)}: acierta {ok}, falla {len(fallos)}")
print(f"fuera de hipotesis (p | M h^-): {fuera_total}, de ellos la formula difiere en {fuera_difiere}")
print(f"sin h^- en la tabla: {sin_h}")
for f in fallos[:20]:
    print("FALLA", f)
