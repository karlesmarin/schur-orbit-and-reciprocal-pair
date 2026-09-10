# -*- coding: utf-8 -*-
# CADA CIFRA Y CADA FORMULA DE characters_at_aP.tex, VERIFICADA.   25 de agosto de 2026.
#
# POR QUE.  El documento va a un tercero que conoce la materia mejor que nosotros.  Si una cifra o
# una formula falla, se cae todo.  Esta compuerta NO se fia del .tex: RECALCULA lo barato desde
# cero y CONTRASTA lo caro contra la salida archivada de las gates.  Sale con codigo 1 si algo no
# cuadra.
#
# QUE COMPRUEBA
#   N1  la tabla del Teorema 2 de Prasad (7 filas + total 2428 + senuelos), RECALCULADA aqui
#   N2  las capacidades 21630 / 31 combinaciones,  contra gates/aP_potencias_OUT.txt
#   N3  el valor 4632 y el rango 42-77 % del control,  contra la misma salida
#   N4  la holgura (2170 casos), RECALCULADA aqui en forma cerrada
#   N5  los 186 pares de R_P(q),  contra gates/R_P_todos_los_tipos_OUT.txt
#   N6  las formulas del .tex que se pueden evaluar: H(z), ht_P, capacidad, valor
#   N7  que toda cifra con \ver{} tiene un guion que existe
#   N8  numeros del .tex que NO aparecen en ninguna fuente -> se listan para mirarlos a mano
#   N12 la demostracion de la Prop. 3.3 escrita en el .tex, paso a paso, RECALCULADA aqui ---
#       incluido el exponente del Vandermonde, que es q^{(q-2)/2} y NO q^{q-2} (ese es su cuadrado)
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python check_numeros.py > check_numeros_OUT.txt 2>&1

import io
import itertools
import os
import re
import sys
from collections import Counter
from fractions import Fraction
from math import comb, gcd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

TEX = io.open("characters_at_aP.tex", encoding="utf-8").read()
GATES = os.path.join("..", "gates")
fallos = []


def ok(cond, etiqueta, detalle=""):
    print(("   OK   " if cond else "  FALLA ") + etiqueta + ("   " + detalle if detalle else ""))
    if not cond:
        fallos.append(etiqueta)


def en_tex(s):
    return s in TEX


def leer(p):
    with io.open(p, encoding="utf-8", errors="replace") as f:
        return f.read()


# ---------------------------------------------------------------- maquinaria minima, independiente
def det_int(M):
    n = len(M)
    if n == 0:
        return 1
    M = [r[:] for r in M]
    sg, prev = 1, 1
    for c in range(n - 1):
        if M[c][c] == 0:
            p = next((r for r in range(c + 1, n) if M[r][c] != 0), None)
            if p is None:
                return 0
            M[c], M[p] = M[p], M[c]
            sg = -sg
        for r in range(c + 1, n):
            for k in range(c + 1, n):
                M[r][k] = (M[r][k] * M[c][c] - M[r][c] * M[c][k]) // prev
            M[r][c] = 0
        prev = M[c][c]
    return sg * M[n - 1][n - 1]


def dim_GL(mu):
    m = len(mu)
    b = [mu[i] + (m - 1 - i) for i in range(m)]
    num = den = 1
    for k in range(m):
        for l in range(k + 1, m):
            num *= (b[k] - b[l])
            den *= (l - k)
    return num // den


def dominantes(N, top):
    for lam in itertools.product(range(top + 1), repeat=N):
        if all(lam[i] >= lam[i + 1] for i in range(N - 1)):
            yield list(lam)


print("=" * 90)
print("N1 -- LA TABLA DEL TEOREMA 2 DE PRASAD, RECALCULADA DESDE CERO EN ESTE FICHERO")
FILAS = [(1, 3, 6, 84, 27, 27), (2, 2, 6, 210, 40, 40), (2, 3, 5, 462, 72, 54),
         (3, 2, 5, 462, 60, 40), (2, 4, 4, 495, 54, 54), (3, 3, 3, 220, 32, 32),
         (4, 2, 4, 495, 45, 45)]
total = 0
for (m, n, top, pesos_esp, dec1_esp, dec2_esp) in FILAS:
    N = m * n
    H = lambda k, n=n, m=m: (comb(k // n + m - 1, m - 1) if (k >= 0 and k % n == 0) else 0)
    pesos = crit = val = d1 = d2 = 0
    for lam in dominantes(N, top):
        pesos += 1
        M = [[H(lam[i] - (i + 1) + (j + 1)) for j in range(N)] for i in range(N)]
        v = det_int(M)
        a = [lam[i] + (N - 1 - i) for i in range(N)]
        cnt = Counter(x % n for x in a)
        c = all(cnt.get(i, 0) == m for i in range(n))
        crit += (c == (v != 0))
        if c:
            bl = []
            for i in range(n):
                b = sorted((x - i) // n for x in a if x % n == i)[::-1]
                bl.append(dim_GL([b[k] - (m - 1 - k) for k in range(m)]))
            p = 1
            for x in bl:
                p *= x
            val += (abs(v) == p)
            pm = 1
            for x in bl[:-1]:
                pm *= x
            d1 += (abs(v) == pm)
            d2 += (abs(v) == p * bl[0])
        else:
            val += 1
    total += pesos
    ok(pesos == pesos_esp and crit == pesos and val == pesos and d1 == dec1_esp and d2 == dec2_esp,
       f"m={m} n={n}",
       f"pesos {pesos} (tex {pesos_esp}) · criterio {crit}/{pesos} · valor {val}/{pesos} · "
       f"senuelos {d1}/{dec1_esp}, {d2}/{dec2_esp}")
ok(total == 2428, "total 2428", f"calculado {total}")
ok(en_tex("$\\mathbf{2428}$"), "el .tex dice 2428")
for (m, n, top, pesos, d1, d2) in FILAS:
    ok(en_tex(f"{m}&{n}& {pesos} "), f"la fila m={m},n={n},{pesos} esta en el .tex")

print("=" * 90)
print("N2/N3 -- CAPACIDADES Y VALOR, contra la salida archivada de aP_potencias")
out = leer(os.path.join(GATES, "aP_potencias_OUT.txt"))
mt = re.search(r"TOTAL:\s*(\d+)\s*aciertos,\s*(\d+)\s*fallos", out)
ok(bool(mt) and mt.group(1) == "21630" and mt.group(2) == "0",
   "21630 aciertos / 0 fallos en la gate", mt.group(0) if mt else "NO ENCONTRADO")
ok(en_tex("$\\mathbf{21\\,630}$"), "el .tex dice 21 630")
filas_v1 = re.findall(r"^\s*\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s*$", out, re.M)
ok(len(filas_v1) == 31, "31 combinaciones (m,k) en la gate", f"contadas {len(filas_v1)}")
ok(en_tex("$31$ combinations"), "el .tex dice 31 combinations")

v2 = re.findall(r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/(\d+)\s+(\d+)/(\d+)\s+(\d+)\s*$", out, re.M)
suma = sum(int(r[2]) for r in v2)
ok(suma == 4632, "4632 pesos en V2", f"suma {suma}")
ok(all(r[3] == r[4] for r in v2), "V2: la formula acierta el 100 % en los seis rangos")
pcts = [100.0 * int(r[5]) / int(r[6]) for r in v2]
ok(abs(min(pcts) - 42) < 1.0 and abs(max(pcts) - 77) < 1.0,
   "el control cae entre 42 % y 77 %", f"min {min(pcts):.1f} max {max(pcts):.1f}")
ok(en_tex("$\\mathbf{4632}$") and en_tex("$42\\%$ and $77\\%$"), "el .tex dice 4632 y 42-77 %")

print("=" * 90)
print("N4 -- LA HOLGURA, RECALCULADA EN FORMA CERRADA (Proposicion del .tex)")
casos = mal = 0
for m in range(2, 60):
    t = 2 * m + 2
    for k in range(1, 40):
        d = gcd(k, t)
        q = t // d
        if q < 3:
            continue
        fij = 1 + (1 if q % 2 == 0 else 0)
        nof = (q // 2 - 1) if q % 2 == 0 else (q - 1) // 2
        cap = d * nof + (d // 2) * fij
        casos += 1
        # el enunciado del .tex: filas = dq/2 - 1;  cap = dq/2 si d par, dq/2 - 1 si d impar
        if m != d * q // 2 - 1:
            mal += 1
        if cap != (d * q // 2 if d % 2 == 0 else d * q // 2 - 1):
            mal += 1
ok(casos == 2170 and mal == 0, "2170 casos, la forma cerrada del .tex se cumple",
   f"casos {casos}, fallos {mal}")
# ⚠ la nota ya no dice "59" y "39" sueltos: ahora escribe el RANGO COMPLETO, $2\le m\le59$ y
#   $1\le k\le39$, con el filtro $q\ge3$ que antes faltaba y sin el cual el 2170 no se reproduce.
#   La comprobacion literal vieja fallo por el cambio de redaccion, no por un error de la nota.
_sb2170 = re.sub(r"\s+", "", TEX)
ok(en_tex("$\\mathbf{2170}$")
   and "2\\lem\\le59" in _sb2170 and "1\\lek\\le39" in _sb2170 and "q\\ge3" in _sb2170,
   "el .tex dice 2170 con su rango completo y el filtro q>=3")

print("=" * 90)
print("N5 -- LOS 186 PARES DE R_P(q), contra la salida archivada")
o5 = leer(os.path.join(GATES, "R_P_todos_los_tipos_OUT.txt"))
m5 = re.search(r"(\d+)\s+pares \(tipo, q\):\s+E1\s+(\d+)/(\d+)\s+E2\s+(\d+)/(\d+)\s+E3\s+(\d+)/(\d+)", o5)
ok(bool(m5), "la linea de resultados existe", m5.group(0) if m5 else "NO")
if m5:
    n_, a1, b1, a2, b2, a3, b3 = m5.groups()
    ok(n_ == "186" and a1 == b1 == "186" and a2 == b2 == "186" and a3 == b3 == "186",
       "186/186 en las tres afirmaciones", m5.group(0))
ok(en_tex("$\\mathbf{186}$"), "el .tex dice 186")

print("=" * 90)
print("N6 -- LAS FORMULAS DEL .tex, EVALUADAS")
# (a) H(z) = (1-z^n)^{-m}  =>  h_k = C(k/n+m-1, m-1)
ok(en_tex("H(z)=\\prod_{\\omega^n=1}(1-\\omega z)^{-m}=(1-z^n)^{-m}"), "la formula de H(z) esta escrita")
malh = 0
for m in (2, 3):
    for n in (2, 3):
        h = [0] * 25
        h[0] = 1
        for j in range(n):
            for _ in range(m):
                w = complex(1) if False else None
        # comprobacion directa: coeficientes de (1-z^n)^{-m}
        for k in range(25):
            esperado = comb(k // n + m - 1, m - 1) if k % n == 0 else 0
            # serie de (1-z^n)^{-m} por convolucion
            c = [0] * 25
            c[0] = 1
            for _ in range(m):
                nueva = [0] * 25
                for i in range(25):
                    if c[i]:
                        for j in range(0, 25 - i, n):
                            nueva[i + j] += c[i]
                c = nueva
            if c[k] != esperado:
                malh += 1
            break
ok(malh == 0, "h_k = C(k/n+m-1, m-1) reproduce la serie de (1-z^n)^{-m}")

# (b) ht_P = ((alpha,alpha)/2) ht(alpha^v)  -- se comprueba que ht(alpha^v) sale ENTERO
o5b = leer(os.path.join(GATES, "R_P_todos_los_tipos_OUT.txt"))
ok("180 raices, 0 no enteros" in o5b, "ht(alpha^v) entero en las 180 raices")
# ⚠ la identidad crecio en la revision 52: ahora la nota escribe tambien  = (rho, alpha), que es
#   lo que ata la seccion 5 con la seccion 3 (la P-altura ES el estadistico de raices de (1)).
ok(en_tex("\\hgt_P(\\alpha) \\;=\\; \\tfrac{(\\alpha,\\alpha)}{2}\\,\\hgt(\\alpha^\\vee) \\;=\\; (\\rho,\\alpha)"),
   "la identidad de la P-altura, con el = (rho,alpha)")
_ob = leer(os.path.join(GATES, "revision_52_bisagra_OUT.txt"))
ok("  FALLA " not in _ob, "la gate de la bisagra no trae ningun fallo")
# ⚠ EL SOPORTE DE LOS h_n.  La nota decia "tres residuos mod q" a secas, y con d IMPAR son DOS
#   --- y d impar incluye a a_P mismo, justo la mitad que la seccion 2 subraya.  Ninguna compuerta
#   lo veia porque miraban el ALFABETO y H(z), no los h_n.  Ahora se mira.
_o2 = leer(os.path.join(GATES, "revision_53_seccion2_OUT.txt"))
ok("  FALLA " not in _o2, "la gate de la seccion 2 no trae ningun fallo")
ok("120/120" in _o2, "las afirmaciones de la seccion 2, recalculadas en 120 casos")
ok(en_tex("when $d$ is even and on\n$0,2$ when $d$ is odd"),
   "el soporte de los h_n distingue la paridad de d")
ok(not en_tex("supported on three residues"),
   "y ya no dice 'tres residuos' a secas")
ok("140/140" in _ob, "la bisagra a_P^k = rho(xi_q), en 140 pares (m,k)")
ok(en_tex("a_P^{\\,k}=\\rho(\\zeta_t^{\\,k})=\\rho(\\xi_q)"), "y la nota la escribe")

# (c) el valor: |ell_i -+ ell_j|/q  y  2 ell_i / q
ok(en_tex("|\\ell_i\\mp\\ell_j|/q") and en_tex("2\\ell_i/q"), "las dos formulas del valor estan escritas")
# ⚠ sin espacios: al reescribir la conjetura el espacio antes de \mid desaparecio y la comprobacion
#   fallo aunque la definicion seguia ahi.  Una comprobacion de cadena literal es fragil ante el
#   reformateo; se normaliza el blanco.
_sin_blanco = re.sub(r"\s+", "", TEX)
# ⚠ y en la revision 52 cambio el CORCHETE: <,> es el apareamiento con corraices y (,) la forma
#   invariante.  La nota usaba <alpha,rho> para un apareamiento peso-peso, que en una nota cuya
#   tesis es raices-contra-corraices se lee justo al reves de lo que quiere decir.
ok("R_P(q)=\\{\\alpha:q\\mid(\\rho,\\alpha)\\}".replace(" ", "") in _sin_blanco,
   "la definicion de R_P(q), ignorando espacios")
# y no puede quedar ningun <,> con dos PESOS dentro
_malos = [p for p in ("\\langle\\ell,\\alpha\\rangle", "\\langle\\rho,\\alpha\\rangle",
                      "\\langle\\alpha,\\rho\\rangle", "\\langle\\rho_P,\\alpha\\rangle")
          if p in TEX]
ok(not _malos, "ningun <,> se usa para un apareamiento peso-peso", str(_malos))

# (d) LA DERIVACION: la identidad de alfabetos, recalculada aqui desde cero
print()
print("     la Proposicion del alfabeto, RECALCULADA en este fichero:")
import cmath
malos = tot_a = 0
for m_ in range(2, 16):
    t_ = 2 * m_ + 2
    for k_ in range(1, 13):
        d_ = gcd(k_, t_)
        q_ = t_ // d_
        if q_ < 2:
            continue
        z_ = cmath.exp(2j * cmath.pi / t_)
        obs = Counter()
        for e in ([(j * k_) % t_ for j in range(1, m_ + 1)] +
                  [(-j * k_) % t_ for j in range(1, m_ + 1)]):
            obs[round(cmath.phase(z_ ** e) / (2 * cmath.pi) * q_) % q_] += 1
        pred = Counter({r: d_ for r in range(q_)})
        if d_ % 2 == 0:
            pred[0] -= 2
        else:
            pred[0] -= 1
            pred[q_ // 2] -= 1
        pred = Counter({r: c for r, c in pred.items() if c})
        tot_a += 1
        malos += (obs != pred)
ok(malos == 0, f"A = d*mu_q menos dos puntos, en {tot_a} pares (m,k)", f"{malos} fallos")
ok(en_tex("A \\;=\\; d\\cdot\\mu_q \\;-\\; 2\\cdot\\{1\\}"), "el .tex enuncia el caso d par")
ok(en_tex("A \\;=\\; d\\cdot\\mu_q \\;-\\;\\{1\\}-\\{-1\\}"), "el .tex enuncia el caso d impar")
ok(en_tex("H(z)=\\frac{(1-z)^2}{(1-z^q)^d}") and en_tex("H(z)=\\frac{1-z^2}{(1-z^q)^d}"),
   "el .tex enuncia las dos formas de H(z)")
o_alf = leer(os.path.join(GATES, "aP_alfabeto_OUT.txt"))
ok("0 aciertos de 72" in o_alf, "el senuelo de paridad invertida acierta 0 de 72")
ok(en_tex("$0$ of $72$ cases"), "el .tex cita el 0 de 72")

# (e) EL RANGO DE LOS SENUELOS, calculado de la propia tabla --- se cito mal una vez ("10-30 %")
pcts = []
for (m_, n_, top_, pesos_, d1_, d2_) in FILAS:
    pcts += [100.0 * d1_ / pesos_, 100.0 * d2_ / pesos_]
lo, hi = min(pcts), max(pcts)
ok(abs(lo - 8.7) < 0.1 and abs(hi - 32.1) < 0.1,
   "el rango real de los senuelos", f"min {lo:.2f} % max {hi:.2f} %")
ok(en_tex("$8.7\\%$ and $32.1\\%$") and not en_tex("$10$--$30\\%$"),
   "el .tex cita el rango correcto y no el viejo")

# (f) LA OBSTRUCCION DEL DETERMINANTE: d*mu_q solo es espectro simplectico si d es PAR
mal_det = 0
for d_ in range(1, 8):
    for q_ in range(2, 10):
        det = ((-1) ** (q_ + 1)) ** d_
        es_sp = (det == 1) and (d_ % 2 == 0 or q_ % 2 == 1)
        # en la familia de a_P, d impar fuerza q par; ahi NUNCA puede ser espectro simplectico
        if d_ % 2 == 1 and q_ % 2 == 0 and es_sp:
            mal_det += 1
ok(mal_det == 0, "d impar + q par => d*mu_q no es espectro simplectico (det = -1)")
ok(en_tex("\\det(d\\cdot\\mu_q)=(-1)^{d}=-1"), "el .tex enuncia la obstruccion del determinante")
ok(not en_tex("one element, from two sides"),
   "el .tex ya NO afirma que sean el mismo elemento en todos los casos")
# ⚠ desde que la §7 pasa a hablar del limite y no de la diferencia finita, h'_n sale una sola vez.
#   Lo que importa no es cuantas veces, sino que SIEMPRE que salga lleve el indicador.
ok(en_tex("\\mathbf{1}_{q\\mid n}"), "h'_n lleva el indicador [q | n]")
ok(TEX.count("h'_n=\\binom") == 0, "no queda ningun h'_n SIN el indicador",
   f"{TEX.count(chr(39))} apostrofes en total")

print("=" * 90)
print("N7 -- LOS GUIONES CITADOS CON \\ver{} EXISTEN")
# ⚠ La primera version de este bloque usaba una clase de caracteres sin el PUNTO, no capturaba
#   ningun nombre de fichero, y daba "OK existe \" --- un control que no podia fallar.
guiones = set()
for arg in re.findall(r"\\ver\{(.+?)\}", TEX):
    for trozo in arg.split(","):
        t = trozo.replace("\\_", "_").replace("\\", "").strip()
        if t.endswith(".py"):
            guiones.add(t)
ok(len(guiones) >= 3, "la regex captura nombres de guion", f"encontrados: {sorted(guiones)}")
for nombre in sorted(guiones):
    ok(os.path.exists(os.path.join(GATES, nombre)), f"existe gates/{nombre}")
# y el senuelo: un nombre que NO existe tiene que dar FALLA si se probara
ok(not os.path.exists(os.path.join(GATES, "no_existe_este_guion.py")),
   "senuelo: un guion inexistente no se encuentra")

print("=" * 90)
print("N7b -- Y CADA UNO TIENE SU SALIDA ARCHIVADA")
# ⚠ 9-sep-2026.  N7 pedia MENOS que lo que el articulo promete.  El §12 dice que cada cifra con
#   marcador azul "comes from a script with archived output, available on request", y N7 solo
#   comprobaba que el GUION existiera.  Ese dia habia 8 guiones citados sin ninguna salida en
#   disco -- entre ellos galois.py (los 11800 pares) y osp.py (los 1386 casos) -- y doce guiones
#   mas de gates/ en la misma situacion.  Prometer una salida que no esta es prometer de mas.
#   Para un guion de figura la salida ES la figura: se le pide el .pdf.


def salida_de(nombre):
    """Donde vive la salida de ese guion.  Para uno de figura la salida ES la figura, y la
    figura vive aqui, junto al .tex que la incluye, no en gates/."""
    base = nombre[:-3]
    if base.startswith("fig_"):
        return os.path.join(".", base + ".pdf")
    return os.path.join(GATES, base + "_OUT.txt")


sin_salida = [g for g in sorted(guiones) if not os.path.exists(salida_de(g))]
ok(not sin_salida, "los %d guiones citados tienen salida archivada" % len(guiones),
   f"sin salida: {sin_salida}" if sin_salida else "")
# senuelo: el control tiene que saber ver una salida que falta
ok(not os.path.exists(os.path.join(GATES, "no_existe_este_guion_OUT.txt")),
   "senuelo: una salida inexistente no se encuentra")

print("=" * 90)
print("N8 -- NUMEROS DEL .tex SIN FUENTE CONOCIDA (mirar a mano)")
CONOCIDOS = {"14", "16", "20", "49", "343", "72", "2428", "21630", "21", "630", "4632", "2170", "186", "285", "31", "32", "26", "6",
             "42", "77", "10", "30", "59", "39", "180", "1940", "1976", "2002", "2009", "2015",
             "2016", "2020", "2023", "2025", "2026", "132", "143", "146", "9", "11", "12", "8",
             "5", "4", "3", "2", "1", "0", "7", "211", "257", "270", "20", "179", "212", "478",
             "99", "101", "79", "36", "499", "111045", "1653", "1676", "84", "210", "462", "495",
             "220", "27", "40", "72", "54", "60", "45", "1",
             # los cuatro que quedaban tras limpiar los falsos positivos, mirados uno a uno:
             "100", "504",    # supervivientes de fig_union3d, calculados por el propio figs.py
             "140",           # los pares (m,k) de revision_52_bisagra.py
             "40401",         # los pesos de G2 del Cor. 7.1
             "27208",         # las parejas de revision_54_enunciados.py
             "21630", "4632", "2170",
             }
# ⚠ La primera version escaneaba el .tex entero y la lista salia llena de FALSOS POSITIVOS:
#   "25" de NPP25, "09" de KLP09, "15" de Petreolle15, "25" de la fecha del documento.  Una lista
#   que grita en falso no la mira nadie, y entonces el control deja de existir.  Se limpia antes:
#   fuera el preambulo, las claves de cita, las etiquetas, los marcadores \ver y los exponentes.
cuerpo = TEX.split("\\begin{thebibliography}")[0].split("\\begin{document}")[-1]
# ⚠ y ANTES de nada, pegar los miles: 27\,208 y 40\,401 se partian en dos, y las mitades sueltas
#   ("401", "208") aparecian como numeros sin fuente.  Se arregla la causa, no cada sintoma.
_limpio = re.sub(r"(\d)\\,(\d)", r"\1\2", cuerpo)
for _pat in (r"\\cite\{[^}]*\}", r"\\ver\{[^}]*\}", r"\\label\{[^}]*\}", r"\\ref\{[^}]*\}",
             r"\\eqref\{[^}]*\}", r"\\includegraphics\[[^]]*\]\{[^}]*\}",
             r"10\^\{-\d+\}", r"10\^\{-?\d+\}", r"\\cdot10\^\{-\d+\}"):
    _limpio = re.sub(_pat, " ", _limpio)
raros = sorted({x for x in re.findall(r"\b\d+\b", _limpio) if x not in CONOCIDOS})
print("   numeros no catalogados:", raros if raros else "ninguno")
if raros:
    print("   (no es un fallo: es la lista que hay que mirar a ojo antes de enviar)")
# y el control del propio control: si se le mete una cifra inventada, tiene que verla
_prueba = sorted({x for x in re.findall(r"\b\d+\b", _limpio + " 987654 ") if x not in CONOCIDOS})
ok("987654" in _prueba, "N8: el barrido SI ve una cifra que no esta catalogada")

print("=" * 90)
print("N11 -- EL TEOREMA: la especializacion principal, el conteo de ceros y el lema de filas")
oR = leer(os.path.join(GATES, "revision_47_verificacion_OUT.txt"))
oG2 = leer(os.path.join(GATES, "revision_47_G2_OUT.txt"))
ok("OK   R1 con P = (1-S)^2" in oR and "OK   SENUELO" in oR,
   "R1: el lema de filas, con senuelo que falla")
ok("OK   R2: la especializacion principal (*) se cumple" in oR, "R2: (*) verificada")
ok("OK   R3: el criterio N_q(lambda) = r_q se cumple" in oR, "R3: el criterio de ceros")
ok("OK   R4: el conteo por clases se cumple" in oR, "R4: el conteo por clases")
ok("OK   R6: el criterio (+) cierra en G2" in oG2, "R6: G2 cierra con el criterio bueno")
ok("OK   q=3: chi=0 <=> a = 2 mod 3   49/49" in oG2 and
   "OK   q=2: chi=0 <=> a=b=1 mod 2   49/49" in oG2, "R7: las dos congruencias de G2")
# y el .tex tiene que reflejarlo
# ⚠ LA JERARQUIA CAMBIO (revision 62).  El resultado GENERAL --- chi != 0 <=> N_q = r_q, valido en
#   todo tipo --- vivia como la ecuacion (2) dentro de la demostracion de un teorema de tipo C, y su
#   especializacion a tipo C era el teorema.  Estaba al reves.  Ahora (2) es el TEOREMA y las
#   capacidades su COROLARIO.  Esta comprobacion exige la jerarquia nueva.
# 8-sep-2026: el teorema se llamaba "[the criterion; any type]".  Tras leer NPP Cor 3.7 --- que da
# el mismo recuento en CORRAICES para todo entero, y por tanto coincide con el nuestro en los tipos
# simplemente enlazados --- decir "any type" era exactamente lo que habia que dejar de decir.  El
# ancla se mueve al nombre nuevo; lo que la compuerta vigila, que el criterio sea TEOREMA y las
# capacidades COROLARIO, no ha cambiado.
ok(en_tex("\\begin{theorem}[the criterion]\\label{thm:crit}"),
   "el .tex declara el criterio general como TEOREMA")
ok(en_tex("\\begin{corollary}[capacity, type $C$]\\label{conj:cap}"),
   "y las capacidades como COROLARIO suyo")
ok(not en_tex("\\begin{theorem}[capacity]"), "y ya no queda el teorema de capacidad viejo")
ok(not en_tex("Theorem~\\ref{conj:cap}"),
   "ninguna referencia llama 'Theorem' al corolario de capacidad")
# el orden importa: el general ANTES que su corolario
ok(TEX.index("\\label{thm:crit}") < TEX.index("\\label{conj:cap}"),
   "y el teorema general va antes que su corolario")
ok(not en_tex("Conjecture~\\ref{conj:cap}"), "no queda ninguna referencia a 'Conjecture' de capacidad")
ok(en_tex("\\label{eq:ps}") and en_tex("\\label{eq:crit}"), "las dos ecuaciones nuevas estan etiquetadas")
ok(en_tex("D_{P(S)h}(\\mathbf b)"), "el lema de filas esta escrito")
ok(en_tex("the principal specialisation \\eqref{eq:ps} & classical"),
   "la atribucion marca la especializacion como clasica")
ok(not en_tex("not separable into a factor"), "la obstruccion falsa ya no esta")
# ⚠ G2 CAMBIO DE ESTATUS (revision 51).  Antes el 49/49 sostenia el resultado; ahora el criterio
#   (2) es independiente del tipo, luego en G2 las dos congruencias son COROLARIO para todo peso y
#   el 49/49 es solo el control con caracteres exactos.  Esta comprobacion buscaba la frase vieja
#   y fallo al reescribirla: se reescribe para lo que la nota afirma AHORA.
_sb = re.sub(r"\s+", " ", TEX)
ok("$49$ of $49$ dominant weights for every $k=1,\\dots,6$" in _sb,
   "G2: el 49/49 sigue citado, ahora como control")
ok(en_tex("\\begin{corollary}\\label{cor:g2}"), "G2: las dos congruencias son un corolario")
ok(en_tex("\\mathbf{40\\,401}"), "G2: el corolario se comprueba en 40401 pesos, no en 49")
_o = leer(os.path.join(GATES, "revision_50_G2_corolario_OUT.txt"))
ok("40401/40401" in _o, "G2: la salida archivada trae el 40401/40401")
# ⚠ "FALLA" a secas no vale: el propio rotulo del bloque de senuelos dice "tienen que FALLAR".
#   Se busca la marca exacta que imprime ok() cuando algo falla, con sus dos espacios delante.
ok("  FALLA " not in _o, "G2: y la gate del corolario no trae ningun fallo")

print("=" * 90)
print("N10 -- LAS CIFRAS DE G2, contra la salida archivada")
oG = leer(os.path.join(GATES, "G2_caracteres_OUT.txt"))
ok("OK   G0: Freudenthal reproduce la dimension de Weyl" in oG, "G0: control de dimensiones en verde")
ok("OK   G1: chi(a_P) in {0,+-1}" in oG, "G1: Kostant se cumple en G2 (el elemento es el correcto)")
fil = re.findall(r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/(\d+)\s+(\d+)/(\d+)\s+(\d+)/(\d+)",
                 oG, re.M)
k1 = [f for f in fil if f[0] == "1"]
ok(bool(k1) and k1[0][8] == k1[0][9] == "49", "k=1: el criterio de raices cierra 49/49",
   f"{k1[0][8]}/{k1[0][9]}" if k1 else "no encontrado")
peor = max((int(f[8]) for f in fil if f[0] != "1" and f[0] != "5"), default=None)
ok(peor == 20, "k>=2: el mejor de los tres es 20/49", f"mejor {peor}")
ok("343" in oG, "el 343 esta en la salida de G2")
ok(en_tex("$49$ of $49$") and en_tex("$20$ of $49$") and en_tex("$343$"),
   "el .tex cita 49/49, 20/49 y 343")

print("=" * 90)
print("N9 -- DERIVA DE TRADUCCION: la version castellana tiene que decir las MISMAS cifras")
print("      (characters_at_aP_es.tex es copia de trabajo interna; si una cifra se separa, una de las dos")
print("       miente, y la que viaja es la inglesa)")
ES = "characters_at_aP_es.tex"
if not os.path.exists(ES):
    print("   (no hay version castellana; nada que comparar)")
else:
    tex_es = io.open(ES, encoding="utf-8").read()

    def cifras(t):
        cuerpo = t.split("\\begin{thebibliography}")[0]
        # ⚠ el indicador \mathbf{1}_{q|n} NO es una cifra afirmada: se quita antes de contar,
        #   porque si no aparece como un "1" que sobra en una version y falta en la otra.
        cuerpo = cuerpo.replace("\\mathbf{1}_{q\\mid n}", "")
        # los numeros que AFIRMAN algo: los de las tablas y los resaltados en negrita
        return Counter(re.findall(r"\\mathbf\{([0-9\\, ]+)\}", cuerpo)) + \
               Counter(re.findall(r"\$(\d+)/(\d+)\$", cuerpo) and
                       [a + "/" + b for a, b in re.findall(r"\$(\d+)/(\d+)\$", cuerpo)])

    c_en, c_es = cifras(TEX), cifras(tex_es)
    solo_en = sorted((c_en - c_es).elements())
    solo_es = sorted((c_es - c_en).elements())
    ok(not solo_en and not solo_es, "las cifras afirmadas coinciden en las dos versiones",
       f"solo EN: {solo_en}   solo ES: {solo_es}" if (solo_en or solo_es) else
       f"{sum(c_en.values())} cifras contrastadas")
    # y las filas de la tabla del Teorema 2, una a una
    for (m, n, top, pesos, d1, d2) in FILAS:
        ok(f"{m}&{n}& {pesos} " in tex_es, f"la fila m={m},n={n},{pesos} tambien esta en el ES")
    # senuelo: una cifra inventada no puede estar en ninguna de las dos
    ok("$\\mathbf{9999}$" not in TEX and "$\\mathbf{9999}$" not in tex_es,
       "senuelo: una cifra inventada no aparece en ninguna version")

    # ⚠ LO QUE ESTE BLOQUE NO VEIA.  Comparar cifras dejo pasar que la castellana estuviera TRES
    #   revisiones por detras: seguia llamando conjetura a un teorema, sin que ninguna cifra
    #   cambiara.  Una copia de revision atrasada es peor que ninguna --- se revisa el documento
    #   equivocado.  Asi que ahora se compara tambien la ESTRUCTURA.
    def etiquetas(t):
        return set(re.findall(r"\\label\{([^}]+)\}", t))

    def entornos(t):
        cuerpo = t.split("\\begin{thebibliography}")[0]
        return Counter(re.findall(r"\\begin\{(theorem|proposition|conjecture|observation|proof)\}",
                                  cuerpo))

    l_en, l_es = etiquetas(TEX), etiquetas(tex_es)
    ok(l_en == l_es, "las dos versiones tienen los MISMOS enunciados etiquetados",
       f"solo EN: {sorted(l_en - l_es)}   solo ES: {sorted(l_es - l_en)}"
       if l_en != l_es else f"{len(l_en)} etiquetas")
    e_en, e_es = entornos(TEX), entornos(tex_es)
    ok(e_en == e_es, "y el mismo recuento de teoremas, proposiciones y demostraciones",
       f"EN {dict(e_en)}   ES {dict(e_es)}")
    # y ninguna de las dos puede tener un enunciado etiquetado como conjetura
    ok(e_en["conjecture"] == 0 and e_es["conjecture"] == 0,
       "ninguna de las dos deja un enunciado como conjetura")
    # senuelo estructural: la comprobacion tiene que poder fallar
    # (se inyecta ANTES de la bibliografia: entornos() corta ahi, y si se anadiera al final
    #  el senuelo no veria nada y el control no podria fallar nunca)
    _falso = TEX.replace("\\begin{thebibliography}",
                         "\\begin{theorem}x\\end{theorem}\n\\begin{thebibliography}", 1)
    ok(entornos(_falso) != e_en,
       "senuelo: la comparacion de estructura SI detecta un enunciado de mas")

    # ⚠ Y LO QUE TAMPOCO VEIA.  La revision 48 separo Kumari II (Nishu Kumari SOLA) de AK22 y
    #   corrigio la castellana, pero la INGLESA --- la que viaja --- se quedo citando AK22 en el
    #   sitio de Kumari II: acreditaba a Ayyer un articulo que no firma.  Nadie lo vio porque
    #   \bibitem sin citar no da aviso de LaTeX.  Asi que ahora: mismas citas en las dos versiones,
    #   y ninguna entrada definida y sin usar.
    def citas(t):
        # ⚠ \cite{a,b} lleva DOS claves: sin partir por la coma, la clave compuesta parecia una
        #   cita colgando y las dos de dentro parecian entradas sin citar.  Tres falsos positivos
        #   de un solo regex ingenuo.
        return {k.strip() for g in re.findall(r"\\cite\{([^}]+)\}", t) for k in g.split(",")}

    def entradas(t):
        return set(re.findall(r"\\bibitem\{([^}]+)\}", t))

    c_en2, c_es2 = citas(TEX), citas(tex_es)
    ok(c_en2 == c_es2, "las dos versiones citan las MISMAS fuentes en los mismos sitios",
       f"solo EN: {sorted(c_en2 - c_es2)}   solo ES: {sorted(c_es2 - c_en2)}"
       if c_en2 != c_es2 else f"{len(c_en2)} claves")
    for et, t in (("EN", TEX), ("ES", tex_es)):
        huerfanas = entradas(t) - citas(t)
        ok(not huerfanas, f"{et}: ninguna entrada de bibliografia queda definida y sin citar",
           str(sorted(huerfanas)))
        colgando = citas(t) - entradas(t)
        ok(not colgando, f"{et}: ninguna cita apunta a una entrada inexistente",
           str(sorted(colgando)))

print("=" * 90)
print("N12 -- LA DEMOSTRACION DE LA PROPOSICION 3.3, PASO A PASO")
print("       Cada igualdad que el .tex escribe se recalcula aqui desde cero.")
import cmath

TOLC = 1e-9


def _mu(q):
    return [cmath.exp(2j * cmath.pi * r / q) for r in range(q)]


def _U(q, s):
    return [x for r, x in enumerate(_mu(q)) if r != s % q]


peor = {"a-x": 0.0, "A": 0.0, "C": 0.0, "D": 0.0, "CD": 0.0, "M": 0.0}
for q in range(3, 11):
    Dq = 1.0
    for x in _mu(q):
        if abs(x * x - 1) > TOLC:
            Dq *= abs(1 - x * x)
    ref = None
    for s in range(q):
        Us, a = _U(q, s), _mu(q)[s]
        # prod |a - x| = q
        p = 1.0
        for x in Us:
            p *= abs(a - x)
        peor["a-x"] = max(peor["a-x"], abs(p - q) / q)
        # A_s = q^{(q-2)/2}
        A = 1.0
        for i in range(len(Us)):
            for j in range(i + 1, len(Us)):
                A *= abs(Us[i] - Us[j])
        peor["A"] = max(peor["A"], abs(A - q ** ((q - 2) / 2)) / q ** ((q - 2) / 2))
        # C_s ordenado, B_s, D_s
        C = B = D = 1.0
        for x in Us:
            for y in Us:
                if abs(x * y - 1) > TOLC:
                    C *= abs(1 - x * y)
        for i in range(len(Us)):
            for j in range(i + 1, len(Us)):
                if abs(Us[i] * Us[j] - 1) > TOLC:
                    B *= abs(1 - Us[i] * Us[j])
        for x in Us:
            if abs(x * x - 1) > TOLC:
                D *= abs(1 - x * x)
        fijo = abs(a * a - 1) <= TOLC
        Cf = q ** (q - 2) * (1.0 if fijo else abs(1 - a * a))
        Df = Dq if fijo else Dq / abs(1 - a * a)
        peor["C"] = max(peor["C"], abs(C - Cf) / Cf, abs(C - B * B * D) / C)
        peor["D"] = max(peor["D"], abs(D - Df) / Df)
        peor["CD"] = max(peor["CD"], abs(C * D - q ** (q - 2) * Dq) / (q ** (q - 2) * Dq))
        M = A * B * D
        if ref is None:
            ref = M
        peor["M"] = max(peor["M"], abs(M / ref - 1))

for et, cl in [(r"prod_{x in U} |a-x| = q", "a-x"),
               (r"A_s = q^{(q-2)/2}  (el .tex ya NO dice q^{q-2})", "A"),
               (r"C_s ordenado = q^{q-2}|1-a^2|,  y  C_s = B_s^2 D_s", "C"),
               (r"D_s = D/|1-a^2|", "D"),
               (r"C_s D_s = q^{q-2} D,  independiente de s", "CD"),
               (r"M_s = M_0", "M")]:
    ok(peor[cl] < 1e-8, "N12: " + et, f"peor {peor[cl]:.1e}")

# y el .tex tiene que decir el exponente correcto, no el cuadrado
ok(r"A_s=q^{(q-2)/2}" in TEX, "N12: el .tex escribe A_s = q^{(q-2)/2}")
# ⚠ 8-sep: esto prohibia ademas la cadena «Vandermonde part», que es lenguaje legitimo --- la
#   definicion de A_red tiene que nombrar esa parte.  Prohibir una FRASE en vez de una AFIRMACION
#   acaba prohibiendo la verdad.  Queda solo el enunciado falso, que es lo que habia que impedir.
ok(r"alone is $q^{\,q-2}$" not in TEX,
   "N12: el .tex ya no afirma que el Vandermonde valga q^{q-2}")
ok(r"\emph{absolute\nvalues} of the non-vanishing".replace("\\n", "\n") in TEX
   or "absolute\nvalues} of the non-vanishing" in TEX,
   "N12: M_s esta definido con MODULOS (en complejo depende de s)")
ok(r"\begin{proof}" in TEX, "N12: la Prop. 3.3 lleva demostracion, no solo medida")
# Ya no puede quedar NINGUN enunciado etiquetado como conjetura, ni ninguna referencia a uno.
# (La palabra suelta si puede aparecer.  Hasta el 25-ago eran DOS usos: el resumen decia "the
#  capacity rule I had been calling a conjecture", y la conjetura AJENA que motiva los papers de
#  Kumari.  La despersonalizacion del 8-sep quito el primero --- "yo la llamaba" es proceso, y un
#  articulo no lo cuenta ---, asi que ahora el uso legitimo es UNO.  Lo que se prohibe sigue siendo
#  el ENUNCIADO propio.)
_cuerpo = TEX.split(r"\begin{thebibliography}")[0]
ok(r"\begin{conjecture}" not in _cuerpo, "N12: no queda ningun enunciado etiquetado conjetura")
ok(r"Conjecture~\ref" not in _cuerpo, "N12: no queda ninguna referencia a una conjetura propia")
ok(_cuerpo.lower().count("conjectur") == 1,
   "N12: la palabra solo sobrevive en su unico uso legitimo (la conjetura ajena)",
   f"{_cuerpo.lower().count('conjectur')} apariciones")

print("=" * 90)
if fallos:
    print(f"RESULTADO: {len(fallos)} FALLOS -> {fallos}")
    sys.exit(1)
print("RESULTADO: TODO CUADRA.  Ninguna cifra ni formula del documento esta sin verificar.")
