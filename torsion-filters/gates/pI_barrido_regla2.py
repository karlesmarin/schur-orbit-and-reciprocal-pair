# -*- coding: utf-8 -*-
# REGLA 2: un instrumento nuevo barre TODO el catalogo el dia en que nace.
# Los instrumentos nacidos el 20 de agosto no se habian pasado por el Paper I.  Aqui se pasan los
# que son GENERICOS --- los que no dependen de la estructura del Paper II.
#
# Lo que se barre, y por que cada uno:
#   A1  ENUNCIADOS SIN FILA DE ATRIBUCION.  Es el control que encontro lem:sdr en el Paper II.
#       Un enunciado sin fila es una reclamacion implicita.
#   A2  ATRIBUCIONES NOMBRADAS SIN CITA.  Es el que encontro Danilov-Grishukhin.  Busca apellidos
#       en el cuerpo que no vayan acompanados de \cite en la misma frase.
#   A3  ETIQUETAS QUE NADA CITA.
#   A4  FRASES DE NARRATIVA DESFASADAS: algo presentado como abierto/medido que el propio paper
#       ya prueba.  Es el que encontro «the one difficulty left in (L1)».
import io
import os
import re
from collections import Counter

P = r"E:\proyectos\Curiosity\research\orbit-pair\paper\orbit_pair.tex"
s = io.open(P, encoding="utf-8").read()
print("Paper I: %s   %d caracteres" % (os.path.basename(P), len(s)))

# ------------------------------------------------------------------ A1
i = s.find("\\section{Attribution}")
if i < 0:
    i = s.find("\\section{What is ours}")
j = s.find("\\section{", i + 10) if i >= 0 else -1
attr = s[i:j] if i >= 0 else ""
print("")
print("=" * 96)
print("A1  ENUNCIADOS SIN FILA EN LA TABLA DE ATRIBUCION")
print("=" * 96)
if not attr:
    print("   *** el Paper I NO TIENE seccion de atribucion --- el control no aplica igual ***")
PAT = re.compile(r"\\begin\{(theorem|proposition|lemma|corollary|conjecture)\}"
                 r"(?:\[([^\]]*)\])?\s*\\label\{([^}]*)\}")
todos = [(m.group(1), (m.group(2) or "").strip(), m.group(3)) for m in PAT.finditer(s)]
sin = [t for t in todos if attr and ("ref{%s}" % t[2]) not in attr]
print("   enunciados: %d   con fila: %d   SIN fila: %d"
      % (len(todos), len(todos) - len(sin), len(sin)))
for tipo, tit, lab in sin[:25]:
    print("      %-12s %-24s %s" % (tipo, lab, tit[:52]))

# ------------------------------------------------------------------ A2
print("")
print("=" * 96)
print("A2  ATRIBUCIONES NOMBRADAS SIN CITA EN LA MISMA FRASE")
print("=" * 96)
APELLIDOS = ["Littlewood", "Macdonald", "Ayyer", "Kumari", "Albion", "Lecouvey", "Koike", "Terada",
             "Ciucu", "Krattenthaler", "Behrend", "Prasad", "Kostant", "Jantzen", "Lusztig",
             "Zolotarev", "Frobenius", "Lerch", "Gauss", "Fuchs", "Schellekens", "Schweigert",
             "Demazure", "Landweber", "Sjamaar", "Harada", "Yacobi", "Sundaram", "Watanabe",
             "Azenhas", "Molev", "Ostrowski", "Zeilberger", "Kassel", "Hidaka", "Itoh",
             "Danilov", "Grishukhin", "Fulkerson", "Gross", "Hall", "Mirsky", "Seymour",
             "Frohmader", "Fulmek", "Jang", "Kwon", "Alexandersson", "Amini", "Pfannerer",
             "Nadimpalli", "Pattanayak", "Polo", "Murty", "Schur", "Pan", "Migotti", "Moree"]
cuerpo = s[:i] if i > 0 else s
frases = re.split(r"(?<=[.;:])\s", cuerpo)
faltan = Counter()
ejemplos = {}
for f in frases:
    if "\\cite" in f:
        continue
    for a in APELLIDOS:
        if re.search(r"\b%s" % a, f):
            faltan[a] += 1
            ejemplos.setdefault(a, " ".join(f.split())[:190])
for a, n in faltan.most_common(20):
    print("   %-16s %3d frases sin \\cite" % (a, n))
    print("        ej: %s" % ejemplos[a])

# ------------------------------------------------------------------ A3
print("")
print("=" * 96)
print("A3  ETIQUETAS DE ENUNCIADO QUE NADA CITA")
print("=" * 96)
n0 = 0
for tipo, tit, lab in todos:
    if len(re.findall(r"\\ref\{%s\}" % re.escape(lab), s)) == 0:
        print("   0 refs : %-12s %-24s %s" % (tipo, lab, tit[:44]))
        n0 += 1
print("   total: %d de %d" % (n0, len(todos)))

# ------------------------------------------------------------------ A4
print("")
print("=" * 96)
print("A4  FRASES QUE PUEDEN HABER QUEDADO DESFASADAS")
print("=" * 96)
SOSPECHA = re.compile(r"(?i)we know of no|know of nothing|is measured, not proved|"
                      r"remains? open|we do not have|the one difficulty|still open|"
                      r"no tool aimed")
for m in SOSPECHA.finditer(s):
    a = max(0, m.start() - 200)
    b = min(len(s), m.end() + 200)
    frag = re.sub(r"\\(ref|eqref|cite)\{[^}]*\}", "", " ".join(s[a:b].split()))
    print("")
    print("   [%s]" % m.group(0))
    print("      ...%s..." % frag[:330])
