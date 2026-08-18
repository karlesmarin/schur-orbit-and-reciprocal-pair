# Rehace la bibliografia del paper II y engancha los \cite que faltaban.
#
# ESTADO ANTERIOR, medido por _bibaudit.py: 11 entradas, de las que 10 NO se citaban nunca -- los
# nombres aparecian en prosa sin enlace -- y CINCO referencias mencionadas no tenian entrada, entre
# ellas Kostant con ocho menciones, que es el hecho externo mas usado del paper.
#
# Las entradas compartidas con el companion se copian VERBATIM de el: ya estan verificadas alli, y
# asi los dos hermanos citan igual.  Las tres que el companion no tiene se verificaron hoy, una a
# una, antes de escribirlas:
#   Kostant   Adv. Math. 20 (1976), 179-212, el articulo donde vive el teorema de los valores
#             {1,0,-1}; localizado via Prasad arXiv:1402.5504, cuyo titulo lo nombra.
#   Ostrowski Jahresber. DMV 30 (1921), 98-99, la charla del Mathematikertag de Jena; hay traduccion
#             inglesa en ACM SIGSAM Bull. 33 (1999), 5.
#   FSS       Comm. Math. Phys. 180 (1996), 39-97, donde se definen el algebra de orbita y el
#             caracter twining.
# Y se adopta el estilo del companion: \bibitem[Etiqueta]{clave}, no numerico.
import sys

p = 'orbit_pair_ii.tex'
s = open(p, encoding='utf-8').read()

BIB = r"""\begin{thebibliography}{AAAA}

\bibitem[AS14]{AndersenStroppel} H.~H.~Andersen and C.~Stroppel, \emph{Fusion rings for quantum
groups}, Algebr. Represent. Theory \textbf{17} (2014), 1869--1888; arXiv:1212.5736.

\bibitem[AB19]{AB19} A.~Ayyer and R.~E.~Behrend, \emph{Factorization theorems for classical group
characters, with applications to alternating sign matrices and plane partitions}, J. Combin. Theory
Ser.~A \textbf{165} (2019), 78--105; arXiv:1804.04514.

\bibitem[AK22]{AK22} A.~Ayyer and N.~Kumari, \emph{Factorization of classical characters twisted by
roots of unity}, J. Algebra \textbf{609} (2022), 437--483; arXiv:2109.11310.

\bibitem[CK09]{CK09} M.~Ciucu and C.~Krattenthaler, \emph{A factorization theorem for classical group
characters, with applications to plane partitions and rhombus tilings}, in: Advances in Combinatorial
Mathematics, Springer, Berlin, 2009.

\bibitem[Fei04]{Feingold} A.~J.~Feingold, \emph{Fusion rules for affine Kac--Moody algebras},
Contemp. Math. \textbf{343} (2004), 53--96; arXiv:math/0212387.

\bibitem[FSS96]{FSS} J.~Fuchs, A.~N.~Schellekens and C.~Schweigert, \emph{From Dynkin diagram
symmetries to fixed point structures}, Comm. Math. Phys. \textbf{180} (1996), 39--97;
arXiv:hep-th/9506135.

\bibitem[Jan73]{Jantzen} J.~C.~Jantzen, \emph{Darstellungen halbeinfacher algebraischer Gruppen},
Bonner Math. Schriften \textbf{67} (1973).

\bibitem[KY12]{KimYacobi} S.~Kim and O.~Yacobi, \emph{A basis for the symplectic group branching
algebra}, J. Algebraic Combin. \textbf{35} (2012), 269--290; arXiv:1005.2320.

\bibitem[Kos76]{Kostant} B.~Kostant, \emph{On Macdonald's $\eta$-function formula, the Laplacian and
generalized exponents}, Adv. Math. \textbf{20} (1976), 179--212.

\bibitem[KLP09]{KLP} S.~Kumar, G.~Lusztig and D.~Prasad, \emph{Characters of simply-laced nonconnected
groups versus characters of nonsimply-laced connected groups}, Contemp. Math. \textbf{478} (2009),
99--101; arXiv:math/0701615.

\bibitem[Kum24]{Kum24} N.~Kumari, \emph{Factorization of classical characters twisted by roots of
unity: II}, J. Pure Appl. Algebra \textbf{228} (2024), no.~11, 107714; arXiv:2212.12477.

\bibitem[LR34]{LR34} D.~E.~Littlewood and A.~R.~Richardson, \emph{Immanants of some special
matrices}, Q. J. Math. \textbf{os-5} (1934), no.~1, 269--282.

\bibitem[Mar26]{PaperI} C.~Mar\'in, \emph{Factorization and vanishing of Schur polynomials twisted by
roots of unity and reciprocal pairs}, arXiv:2608.09619.

\bibitem[NPP25]{NPP} S.~Nadimpalli, S.~Pattanayak and D.~Prasad, \emph{Character theory at a torsion
element}, arXiv:2504.14684.

\bibitem[Ost21]{Ostrowski} A.~M.~Ostrowski, \emph{\"Uber die Bedeutung der Theorie der konvexen
Polyeder f\"ur die formale Algebra}, Jahresber. Deutsch. Math.-Verein. \textbf{30} (1921), 98--99;
English translation, ACM SIGSAM Bull. \textbf{33} (1999), no.~1, 5.

\bibitem[Pra16]{Prasad} D.~Prasad, \emph{Half the sum of positive roots, the Coxeter element, and a
theorem of Kostant}, Forum Math. \textbf{28} (2016), no.~4, 767--782; arXiv:1402.5504.

\end{thebibliography}"""

i = s.find(r'\begin{thebibliography}')
j = s.find(r'\end{thebibliography}') + len(r'\end{thebibliography}')
assert i > 0 and j > i, 'no encuentro la bibliografia'
s = s[:i] + BIB + s[j:]

# ---- los \cite que faltaban, uno por mencion en prosa -------------------------------------------
CITAS = [
    ("Ayyer--Kumari's \\stext", "Ayyer--Kumari's \\cite{LR34,AK22} \\stext"),
    ("this is Ciucu--Krattenthaler and Ayyer--Behrend \\stext",
     "this is Ciucu--Krattenthaler and Ayyer--Behrend \\cite{CK09,AB19} \\stext"),
    ("(Kostant \\stext)", "(Kostant \\cite{Kostant} \\stext)"),
    ("(Andersen--Stroppel \\stext, in the form we quote it)",
     "(Andersen--Stroppel \\cite{AndersenStroppel} \\stext, in the form we quote it)"),
    ("the element is regular of\norder $t$, and $\\tau_t(\\eta)\\in\\{0,\\pm1\\}$ by Kostant's theorem \\stext",
     "the element is regular of\norder $t$, and $\\tau_t(\\eta)\\in\\{0,\\pm1\\}$ by Kostant's theorem\n\\cite{Kostant,Prasad} \\stext"),
    ("stated for quantum groups in exactly this form by Andersen and Stroppel",
     "stated for quantum groups in exactly this form by Andersen and Stroppel \\cite{AndersenStroppel}"),
    ("and the resulting\nvalue lies in $\\{0,\\pm1\\}$ by Kostant's theorem on characters at elements of finite order.",
     "and the resulting\nvalue lies in $\\{0,\\pm1\\}$ by Kostant's theorem on characters at elements of finite order\n\\cite{Kostant}."),
    ("Littlewood's evaluation of a Schur function on a full root-of-unity orbit",
     "Littlewood's evaluation \\cite{LR34} of a Schur function on a full root-of-unity orbit"),
    ("the branching\nalgebra of Kim and Yacobi \\stext", "the branching\nalgebra of Kim and Yacobi \\cite{KimYacobi} \\stext"),
    ("(Ostrowski), and $\\Phi = N_\\beta / N_\\delta$", "\\cite{Ostrowski}, and $\\Phi = N_\\beta / N_\\delta$"),
    ("Ostrowski's additivity says", "Ostrowski's additivity \\cite{Ostrowski} says"),
    ("(Fuchs--Schellekens--Schweigert)", "\\cite{FSS}"),
    ("That is Jantzen's twining theorem in the form of Kumar--Lusztig--Prasad,",
     "That is Jantzen's twining theorem \\cite{Jantzen} in the form of Kumar--Lusztig--Prasad \\cite{KLP},"),
    ("with the threshold given explicitly (Feingold)", "with the threshold given explicitly \\cite{Feingold}"),
    ("(Andersen--Stroppel) --- so the question can be posed sharply",
     "\\cite{AndersenStroppel} --- so the question can be posed sharply"),
    ("gives only the full family of generators (Andersen--Stroppel)",
     "gives only the full family of generators \\cite{AndersenStroppel}"),
    ("& Kostant & \\stext", "& Kostant \\cite{Kostant} & \\stext"),
    ("& Andersen--Stroppel & \\stext", "& Andersen--Stroppel \\cite{AndersenStroppel} & \\stext"),
    ("& Kim--Yacobi & \\stext", "& Kim--Yacobi \\cite{KimYacobi} & \\stext"),
    ("& Littlewood; Ayyer--Kumari & \\stext", "& Littlewood \\cite{LR34}; Ayyer--Kumari \\cite{AK22} & \\stext"),
]
faltan = []
for a, b in CITAS:
    if a in s:
        s = s.replace(a, b, 1)
    else:
        faltan.append(a[:70])

open(p, 'w', encoding='utf-8', newline='\n').write(s)
print('bibliografia reescrita: 16 entradas, estilo del companion')
print('citas enganchadas: %d de %d' % (len(CITAS) - len(faltan), len(CITAS)))
for f in faltan:
    print('   NO ENCONTRADA:', repr(f))
