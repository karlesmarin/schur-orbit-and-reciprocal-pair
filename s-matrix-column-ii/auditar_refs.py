# -*- coding: utf-8 -*-
r"""auditar_refs.py -- la auditoria de referencias que un arbitro hace en la primera lectura.

Tres comprobaciones, y las tres han cazado algo alguna vez:

 (1) todo \cite tiene su \bibitem, y todo \bibitem se cita;
 (2) TODO APELLIDO QUE APARECE EN LA PROSA tiene entrada en la bibliografia, o esta en la
     lista de eponimos estandar (Newton, Hasse-Weil, ...) que no la necesitan.  Esta es la
     trampa registrada: un nombre en la primera frase sin \bibitem;
 (3) las fechas sueltas en la prosa -- '1786', '1801' -- que invitan a una cita y a menudo
     no la tienen.

    python auditar_refs.py

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import io
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = os.path.dirname(os.path.abspath(__file__))
TEX = io.open(os.path.join(AQUI, "P3b.tex"), encoding="utf-8").read()

# epónimos que nombran un resultado estandar y no piden entrada propia
EPONIMOS = {
    "Newton", "Hasse", "Weil", "Jacobi", "Legendre", "Frobenius", "Galois",
    "Gorenstein", "Weierstrass", "Riemann", "Euler", "Taylor", "Cohen", "Macaulay",
    "Wieferich", "Fermat", "Bring", "Jerrard", "Sato", "Tate", "Clebsch", "Segre",
    "Cayley", "Kneser", "Petersen", "Dirichlet", "Vieweg", "Braunschweig", "Leipzig",
}

cuerpo = TEX.split(r"\begin{thebibliography}")[0]
biblio = TEX[TEX.find(r"\begin{thebibliography}"):]

cites = set(re.findall(r"\\cite\{([^}]+)\}", TEX))
bibs = set(re.findall(r"\\bibitem\{([^}]+)\}", TEX))
print("(1) CLAVES")
print("    citados sin bibitem :", sorted(cites - bibs) or "ninguno")
print("    bibitems sin citar  :", sorted(bibs - cites) or "ninguno")

# apellidos de la bibliografia, para saber quien SI esta
en_biblio = set(re.findall(r"[A-Z][a-z\u00e0-\u00ff]{2,}", re.sub(r"\\emph\{[^}]*\}", "", biblio)))

# apellidos en la prosa: palabras capitalizadas fuera de matematicas, comandos y comienzos
prosa = re.sub(r"\$[^$]*\$", " ", cuerpo)
prosa = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?(\{[^{}]*\})?", " ", prosa)
prosa = re.sub(r"%.*", " ", prosa)
cand = {}
for m in re.finditer(r"(?<![.!?]\s)(?<!^)\b([A-Z][a-z\u00e0-\u00ff]{3,})\b", prosa, re.M):
    cand.setdefault(m.group(1), 0)
    cand[m.group(1)] += 1

PALABRAS = {"Part", "Section", "Theorem", "Proposition", "Corollary", "Their", "These",
            "This", "That", "What", "Where", "When", "Here", "Every", "Each", "From",
            "With", "Since", "Fix", "Let", "The", "Left", "Right", "Both", "Away",
            "Reading", "Expanding", "Counting", "Pulling", "Dividing", "Comparing",
            "Applying", "Writing", "Finally", "Statements", "Infinite", "Outside",
            "Numerically", "Deciding", "Between", "Measured", "January", "Which",
            "Stopping", "Passing", "Additive", "Extending", "Away", "Below", "Above"}
sospechosos = {k: v for k, v in cand.items()
               if k not in PALABRAS and k not in EPONIMOS and k not in en_biblio}
print("")
print("(2) APELLIDOS EN LA PROSA sin entrada en la bibliografia ni en la lista de epónimos")
if sospechosos:
    for k, v in sorted(sospechosos.items(), key=lambda x: -x[1]):
        print("    %-18s x%d" % (k, v))
else:
    print("    ninguno")

print("")
print("(5) FIGURAS Y ECUACIONES HUERFANAS: una figura que ninguna frase nombra flota sola")
for tipo, pat in (("figura", r"fig:[A-Za-z]+"), ("ecuacion", r"eq:[A-Za-z]+"),
                  ("teorema", r"(?:thm|prop|cor):[A-Za-z0-9]+")):
    etiquetas = set(re.findall(r"\\label\{(" + pat + r")\}", TEX))
    for e in sorted(etiquetas):
        refs = len(re.findall(r"\\(?:eq)?ref\{" + re.escape(e) + r"\}", TEX))
        marca = "ok" if refs else "*** HUERFANA ***"
        if tipo == "figura" or not refs:
            print("    %-9s %-14s referida %d vez/veces   %s" % (tipo, e, refs, marca))

print("")
print("(4) EL HILO DE HISTORIA: cada ancla, su fecha, y la SECCION que justifica.")
print("    La regla de la casa: el parrafo historico que no nombre ninguna seccion, fuera.")
HILO = [
    (1786, "Bring", "sec:bring", "la forma normal del estrato extremo"),
    (1801, "Gauss", "sec:centro", "el semisistema, y su lema da el (2|p) del recuento"),
    (1894, "Dedekind", "sec:intro", "el conductor: el invariante que mide la identidad"),
    (1965, "Rosa--Znam", "sec:centro", "el recuento ya existia: la retirada de novedad"),
    (1970, "Kunz", "sec:centro", "simetrico <=> Gorenstein, sin monomialidad"),
    (1973, "Sloane", "sec:centro", "la sucesion catalogada y sin formula"),
    (1988, "Verlinde", "sec:intro", "la columna como dimensiones cuanticas"),
    (2021, "Korchmaros et al.", "sec:bring", "la curva de Bring generalizada y su genero"),
]
for ano, quien, sec, papel in HILO:
    marca = r"\label{%s}" % sec
    print("    %-5d %-18s -> %-14s %s   %s"
          % (ano, quien, sec, papel, "ok" if marca in TEX else "SECCION NO EXISTE"))
huecos = [s for s in ("sec:defecto", "sec:ley") if
          not any(x[2] == s for x in HILO)]
print("    secciones SIN ancla historica: %s" % (huecos or "ninguna"))
print("    (son las dos tecnicas: de donde sale el defecto, y la ley.  No se les fuerza una.)")

print("")
print("(3) FECHAS SUELTAS en la prosa, y si la frase lleva cita cerca")
for m in re.finditer(r"\b(1[6-9]\d\d|20[0-2]\d)\b", cuerpo):
    a, b = max(0, m.start() - 170), min(len(cuerpo), m.end() + 90)
    frag = cuerpo[a:b]
    tiene = bool(re.search(r"\\cite\{", frag))
    print("    %-6s %-4s ...%s..." % (m.group(1), "cita" if tiene else "SIN",
                                      " ".join(frag.split())[-105:]))
