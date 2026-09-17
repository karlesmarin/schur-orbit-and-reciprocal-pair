# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# Reune el indice [O+:A_m(q)] de todas las casillas calculadas (logs de los barridos, en este directorio) y lo
# descompone por primos, con la familia (h, h+1, h+2) de cada primo del soporte.  Escribe
# datos_indice.json, que leen las figuras.  Casillas: fila del reticulo (idx), "conductor trivial"
# (idx = 1) o H2|/H2NADA| (a por primo).  Solo casillas con resultado; nada se interpola.
import glob
import json
import os
import re
import sys
from math import gcd

GATES = os.path.dirname(os.path.abspath(__file__))
PFX = re.compile(r"^(?:(?:sage:|\.\.\.\.:)\s*)*")
FILA = re.compile(r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(SI|no)\s+(.*)$")
CAB = re.compile(r"--- casilla (\S+?)\s+---")


def primos(n):
    out, d = [], 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def vp(x, p):
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v


def familia(m, q, p):
    v = vp(q, p)
    n = q // p ** v
    h = 2 * m
    for nombre, x in (("h", h), ("h+1", h + 1), ("h+2", h + 2)):
        if n >= 3 and x % n == 0:
            return nombre, n
    return None, n


def main():
    idx = {}
    # celdas_m*: filas del reticulo y conductor trivial
    for ruta in sorted(glob.glob(os.path.join(GATES, "celdas_m*_OUT.txt"))):
        t = open(ruta, encoding="utf-8", errors="replace").read()
        cabs = list(CAB.finditer(t))
        for k, mt in enumerate(cabs):
            b = t[mt.end(): cabs[k + 1].start() if k + 1 < len(cabs) else len(t)]
            mq = mt.group(1).split(":")
            m, q = int(mq[0]), int(mq[1])
            got = None
            for ln in b.splitlines():
                f = FILA.match(PFX.sub("", ln).rstrip())
                if f and int(f.group(1)) == m and int(f.group(2)) == q:
                    got = int(f.group(4))
                    break
            if got is None and "conductor trivial" in b:
                got = 1
            if got is not None:
                idx[(m, q)] = got
    # h2 / r2: a por primo
    for nombre in ("h2_OUT.txt", "r2_OUT.txt"):
        ruta = os.path.join(GATES, nombre)
        if not os.path.exists(ruta):
            continue
        partes = {}
        for ln in open(ruta, encoding="utf-8", errors="replace"):
            ln = PFX.sub("", ln.strip())
            if ln.startswith("H2|"):
                _, m, q, p, a, b, dm, forma, c1, c2 = ln.split("|")
                partes.setdefault((int(m), int(q)), 1)
                partes[(int(m), int(q))] *= int(p) ** int(a)
            elif ln.startswith("H2NADA|"):
                _, m, q = ln.split("|")[:3]
                partes[(int(m), int(q))] = 1
        for k, v in partes.items():
            idx.setdefault(k, v)
    celdas = []
    fam_mixta = 0
    for (m, q), I in sorted(idx.items()):
        if not (2 * m < q):
            continue
        porp = []
        for p in primos(q):
            a = vp(I, p)
            fam, n = familia(m, q, p)
            porp.append({"p": p, "a": a, "familia": fam, "qprima": n})
        resto = I
        for p in primos(q):
            while resto % p == 0:
                resto //= p
        fams = set(d["familia"] for d in porp if d["a"] > 0)
        if len(fams) > 1:
            fam_mixta += 1
        celdas.append({"m": m, "q": q, "indice": I, "resto_fuera_de_q": resto, "primos": porp})
    az = [{"m": m, "q": q} for m in range(1, 11) for q in (2 * m + 1, 2 * m + 2)]
    out = {"celdas": celdas, "A_igual_Z": az}
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "datos_indice.json"), "w"), indent=1)
    malos = [c for c in celdas if c["resto_fuera_de_q"] != 1]
    print("casillas: %d   con indice > 1: %d   familias mixtas: %d   primos fuera de q (tiene que ser 0): %d"
          % (len(celdas), sum(1 for c in celdas if c["indice"] > 1), fam_mixta, len(malos)))
    sop_mal = 0
    for c in celdas:
        for d in c["primos"]:
            if d["a"] > 0 and d["familia"] is None and d["qprima"] not in (1, 2, 3, 4, 6):
                sop_mal += 1
    print("primos en el soporte fuera de h, h+1, h+2 con q' no en {1,2,3,4,6} (tiene que ser 0): %d" % sop_mal)
    print("m presentes:", sorted(set(c["m"] for c in celdas)), "  q max:", max(c["q"] for c in celdas))


if __name__ == "__main__":
    main()
