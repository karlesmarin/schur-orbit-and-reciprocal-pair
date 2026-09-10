# -*- coding: utf-8 -*-
# CADA AFIRMACION SOBRE UN ARTICULO AJENO, CONTRA SU FUENTE EN DISCO.   25 de agosto de 2026.
#
# POR QUE.  El documento le atribuye cosas concretas a seis articulos, cuatro de ellos suyos.  Si
# una atribucion esta mal, el destinatario lo ve en el primer parrafo.  Esta compuerta no acepta
# "lo lei": exige que la frase este EN EL FICHERO.
#
# ⚠ AVISO DE SITIO.  Lo que queda fuera de los paquetes que viajan queda fuera POR SITIO y no por
#    lista: la regla y su motivo estan en _privado/README.md.  El control C3 de abajo comprueba
#    que este directorio no se cuele en ninguno de ellos.
#
# QUE COMPRUEBA
#   C1  las frases literales que el .tex atribuye, contra los .txt/.pdf de _papers/
#   C2  los SENUELOS: frases plausibles que NO estan en las fuentes deben salir ausentes
#   C3  que este directorio no esta dentro de ningun paquete que viaje
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python check_citas.py > check_citas_OUT.txt 2>&1

import io
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
TEX = io.open("characters_at_aP.tex", encoding="utf-8").read()
PAP = os.path.join("..", "_papers")
fallos = []


def norm(s):
    s = s.replace("’", "'").replace("“", '"').replace("”", '"')
    s = s.replace("ﬁ", "fi").replace("ﬂ", "fl").replace("−", "-")
    return re.sub(r"\s+", " ", s).strip().lower()


def fuente(nombre):
    p = os.path.join(PAP, nombre)
    if p.endswith(".pdf"):
        import fitz
        d = fitz.open(p)
        return norm("".join(d[i].get_text() for i in range(d.page_count)))
    return norm(io.open(p, encoding="utf-8", errors="replace").read())


FUENTES = {
    "prasad2016": "_prasad_israe.txt",
    "klp": "_klp.txt",
    "polo": "polo_2504.09204.txt",
    "lpa": "lubeck_prasad_ayyer_1912.08576.pdf",
    # ⚠ FALTABA.  La nota afirmaba en el pie de la figura 1 que Ayyer-Kumari no mencionan a Kostant
    #   "ni una vez", y su articulo no estaba ni siquiera en esta lista: ninguna compuerta podia
    #   comprobarlo.  Su introduccion dice "motivated by a celebrated result of Kostant [Kos76]".
    #   Una fuente sobre la que la nota afirma algo TIENE que estar aqui.
    "ak": "arxiv_2109.11310.txt",
    "npp": "NPP25.txt",
    # ⚠ FALTABAN, y es el mismo fallo de la vez anterior con Ayyer-Kumari.  El 9-sep-2026 el
    #   keybox de la §7 y la §13 pasaron a afirmar cosas LITERALES sobre estos dos articulos --
    #   que la clase 0 lleva el factor de raiz larga y es de tipo C, que con l par aparece un
    #   Vandermonde de tipo B_{r_p}, y que su Prop. 3.2.5 dice que los coeficientes no son de
    #   ramificacion.  Ninguna estaba en esta lista, asi que esta compuerta pasaba en verde SIN
    #   MIRARLAS.  Una compuerta que no lleva la afirmacion en su lista no la comprueba.
    "lec_a": "lecouvey_math_0607038.txt",
    "lec_b": "lecouvey_math_0703514.txt",
}
TXT = {}
for k, v in FUENTES.items():
    try:
        TXT[k] = fuente(v)
        print(f"   fuente {k}: {v}  ({len(TXT[k])} caracteres)")
    except Exception as e:
        print(f"  FALLA  no se pudo leer {v}: {e}")
        fallos.append(v)

print("=" * 92)
print("C1 -- LAS FRASES QUE EL DOCUMENTO ATRIBUYE, CONTRA SU FUENTE")

CITAS = [
    ("prasad2016", "is not identically zero if and only if the mn-tuple of integers",
     "Thm 2: el criterio, literal"),
    ("prasad2016", "represents each residue class in z/n exactly m-times",
     "Thm 2: 'exactly m-times'"),
    ("prasad2016", "generalizing a theorem of kostant on the character values at the coxeter element",
     "el resumen de 2016: 'generalizing a theorem of Kostant'"),
    ("prasad2016", "takes one of the values 1, 0, -1",
     "su Teorema 1 = Kostant: valores {1,0,-1}"),
    ("klp", "tr(tσ, ˆπ) = tr(ϖ(t), πσ)".lower().replace(" ", " "),
     "KLP Thm 1.1 (forma con simbolos; si falla se prueba la variante)"),
    ("klp", "characters of simplylaced nonconnected", "KLP: el titulo"),
    # la fuente dice "be the L-dual", no "is the L-dual" --- corregido tras fallar la compuerta
    ("klp", "be the l-dual of the ﬁxed point sub-".replace("ﬁ", "fi"),
     "KLP: G_sigma es el dual del subgrupo de puntos fijos"),
    ("klp", "h = 0 for all η but one and h = 1 for one η".lower(),
     "KLP §2: h=1 para exactamente una orbita en A_2n"),
    ("polo", "the set of roots whose height is a multiple of", "Polo: R(m) por ALTURA"),
    ("polo", "by weyl dimension formula, dm is the dimension of the", "Polo: d_m es una dimension de Weyl"),
    ("polo", "principal elements of order", "Polo: elementos PRINCIPALES"),
    ("lpa", "the dimension of an irreducible representation of bn", "LPA Thm 6.1: +- una dimension"),
    ("lpa", "with an appendix by arvind ayyer", "LPA: Ayyer firma el APENDICE, no el articulo"),
    ("lpa", "equations", "LPA: cita a Littlewood por ecuaciones"),
    ("lpa", "7.3.1 and 7.3.2 of page 132", "LPA: Littlewood ec. 7.3.1-7.3.2, pagina 132"),
    ("lpa", "pages 143-146", "LPA: la forma para S_n, paginas 143-146"),
    ("lpa", "a character relationship on gln(c); israel journal, vol. 211, (2016), 257-270",
     "LPA: [P] ES el Israel J. Math. 211 (2016) 257-270"),
    # --- Lecouvey.  Lo que el keybox de la §7 y la §13 le atribuyen desde el 9-sep-2026.
    ("lec_a", "parabolic kazhdan-lusztig polynomials", "Lec07a: el titulo"),
    ("lec_a", "for g = sp2n", "Lec07a: la seccion 3.2.2 ES de Sp_2n, que es donde miramos"),
    ("lec_a", "(1 -xrxs)",
     "Lec07a: la clase 0 lleva el factor de raiz LARGA -- por eso es tipo C, no tipo A"),
    ("lec_a", "the half sum of positive roots is equal to",
     "Lec07a: la semisuma que identifica el bloque de tipo B"),
    ("lec_a", "in type brp", "Lec07a: ...y la nombra, tipo B_{r_p}"),
    ("lec_a", "these roots do not belong to the root lattice associated",
     "Lec07a: esas raices NO estan en el reticulo de Sp_2n -- el obstaculo, dicho por el"),
    ("lec_a", "there cannot exist an analogue of theorem 3.2.3 when",
     "Lec07a: con l par NO hay analogo de su Thm 3.2.3"),
    ("lec_a", "suppose g = sp2n and", "Lec07a: su Prop. 3.2.5 es de Sp_2n"),
    ("lec_a", "cannot be interpreted as branching coe",
     "Lec07a Prop. 3.2.5: los coeficientes no son de ramificacion"),
    ("lec_a", "and have signs alternatively positive and negative",
     "Lec07a Prop. 3.2.5: ...y sus signos alternan"),
    ("lec_a", "consider a partition µ of length n and ℓ= 2p −1",
     "Lec07a Thm 3.2.3: es el caso l IMPAR, que es el que cierra como identidad"),
    ("lec_b", "stabilized plethysms for the classical lie groups", "Lec07b: el titulo"),
    ("lec_b", "let g be a symplectic or orthogonal",
     "Lec07b Thm 4.5.1: vale para simplectico u ortogonal"),
    ("lec_b", "corresponding to the restriction to certain levi subgroups",
     "Lec07b: en el rango estable SI son coeficientes de ramificacion a un Levi"),
    # --- Ayyer-Kumari Thm 2.11.  ⚠ 9-sep-2026: el keybox de la §7 llego a decir que la
    #     deleccion de las dos letras era "lo que compra el caso par", frente a la Prop. 3.2.5 de
    #     Lecouvey.  ERA FALSO y mezclaba dos objetos: su proposicion habla de los COEFICIENTES de
    #     una descomposicion, y el VALOR si factoriza en orden par SIN quitar ninguna letra --
    #     esto es lo que lo demuestra.  El keybox lo dice ahora, y esta compuerta lo respalda.
    ("ak", "indexing an irreducible repre- sentation of sp2tn",
     "AK22 Thm 2.11: es de Sp_{2tn}, con el alfabeto COMPLETO"),
    ("ak", "if coret(λ) is a symplectic t-core with rank r",
     "AK22 Thm 2.11(2): la rama que SI factoriza"),
    ("ak", "so λ( t 2 −1)(xt)",
     "AK22 Thm 2.11(2): y su factor ORTOGONAL IMPAR, el del caso t par"),
]
for clave, frase, etiqueta in CITAS:
    if clave not in TXT:
        continue
    hay = norm(frase) in TXT[clave]
    print(("   OK   " if hay else "  FALLA ") + etiqueta)
    if not hay:
        fallos.append(etiqueta)

print("=" * 92)
print("C1b -- LAS AFIRMACIONES NEGATIVAS DE LA NOTA, CONTRA LAS FUENTES")
print("      ⚠ Aqui es donde esta compuerta fallo.  Probaba los senuelos que YO recorde, no las")
print("        afirmaciones que la NOTA hace.  El pie de la figura 1 decia que Ayyer-Kumari no")
print("        mencionan a Kostant 'ni una vez', y su introduccion dice literalmente 'motivated")
print("        by a celebrated result of Kostant [Kos76]'.  Una afirmacion negativa sobre el")
print("        articulo de un tercero es la peor de equivocar, y nadie la estaba mirando.")
TEX_NOTA = ""
try:
    TEX_NOTA = io.open("characters_at_aP.tex", encoding="utf-8").read()
except Exception:
    pass

# Toda afirmacion de la forma "X never mentions Y" que la nota haga tiene que tener aqui su prueba.
# Si se anade una nueva a la nota y no se anade aqui, el bloque C1c de abajo lo detecta.
NEGATIVAS = [
    ("ak", "kostant", False,
     "Ayyer-Kumari SI mencionan a Kostant (por eso la nota ya no dice lo contrario)"),
    ("ak", "coxeter", None, "informativo: ¿mencionan el elemento de Coxeter?"),
    ("ak", "torsion", None, "informativo: ¿mencionan elementos de torsion?"),
]
for clave, palabra, esperado_ausente, etiqueta in NEGATIVAS:
    if clave not in TXT:
        print(f"   (sin fuente para {clave}: no se puede comprobar '{palabra}')")
        continue
    hay = norm(palabra) in TXT[clave]
    if esperado_ausente is None:
        print(f"   INFO   {etiqueta}: {'SI aparece' if hay else 'no aparece'}")
        continue
    cond = (not hay) if esperado_ausente else hay
    print(("   OK   " if cond else "  FALLA ") + etiqueta)
    if not cond:
        fallos.append(etiqueta)

# C1c: la nota no puede contener ninguna afirmacion negativa que no este probada arriba
for patron, etq in [("never mention", "la nota no afirma 'never mention' sin prueba en C1b"),
                    ("not once", "la nota no afirma 'not once' sin prueba en C1b"),
                    ("nobody on that line", "la nota no afirma 'nobody on that line'")]:
    hay = patron in TEX_NOTA.lower()
    print(("   OK   " if not hay else "  FALLA ") + etq)
    if hay:
        fallos.append(etq)

print("=" * 92)
print("C2 -- SENUELOS: frases plausibles que NO deben estar (si aparecen, el test no discrimina)")
SENUELOS = [
    ("prasad2016", "the second distinguished element a_p", "Prasad 2016 NO habla de a_P"),
    ("prasad2016", "folded classes", "Prasad 2016 NO habla de clases plegadas"),
    ("polo", "kostant's second element", "Polo NO menciona el segundo elemento de Kostant"),
    ("polo", "each residue class exactly", "Polo NO usa el criterio de residuos de 2016"),
    ("klp", "symplectic character at a_p", "KLP NO menciona a_P"),
    ("lpa", "powers of a_p", "LPA NO habla de potencias de a_P"),
    # --- Lecouvey.  Los senuelos aqui protegen la FRONTERA de lo que le atribuimos: el keybox
    #     dice que la huella de bloques es suya PARA EL ALFABETO COMPLETO, y que lo nuestro es
    #     que sobreviva al borrado de dos letras y en un punto.  Si alguna de estas apareciera,
    #     esa frontera estaria mal puesta y le estariamos atribuyendo de menos.
    ("lec_a", "deleted alphabet", "Lecouvey NO trata el alfabeto BORRADO -- ese es el nuestro"),
    ("lec_a", "short of two letters", "Lecouvey NO quita dos letras del alfabeto"),
    ("lec_a", "kostant's second element", "Lecouvey NO menciona el segundo elemento de Kostant"),
    ("lec_a", "a_p", "Lecouvey NO habla de a_P"),
    ("lec_a", "the value of the character at", "Lecouvey descompone, no evalua en un punto"),
    ("lec_b", "deleted alphabet", "Lec07b tampoco borra letras"),
    ("lec_b", "cannot be interpreted as branching coe",
     "la imposibilidad con l par esta en Lec07a, NO en Lec07b -- no confundir de articulo"),
]
for clave, frase, etiqueta in SENUELOS:
    if clave not in TXT:
        continue
    hay = norm(frase) in TXT[clave]
    print(("   OK   " if not hay else "  FALLA ") + etiqueta + ("" if not hay else "  <-- APARECE"))
    if hay:
        fallos.append("senuelo " + etiqueta)

print("=" * 92)
print("C3 -- SITIO: este directorio no viaja")
raiz = os.path.join("..", "..", "..")
publicados = []
for pk in (os.path.join("..", "paper", "anc"), os.path.join("..", "paper2", "anc")):
    if os.path.isdir(pk):
        for r, _, fs in os.walk(pk):
            for f in fs:
                publicados.append(os.path.join(r, f))
malos = [p for p in publicados if "note_prasad2" in p.replace("\\", "/")]
print(f"   ficheros en los dos paquetes: {len(publicados)}")
print(("   OK   " if not malos else "  FALLA ") + "note_prasad2/ no esta en ningun paquete"
      + ("" if not malos else f"  {malos[:3]}"))
if malos:
    fallos.append("note_prasad2 en un paquete publicado")

print("=" * 92)
if fallos:
    print(f"RESULTADO: {len(fallos)} FALLOS -> {fallos}")
    sys.exit(1)
print("RESULTADO: todas las atribuciones estan en su fuente, y los senuelos no.")
