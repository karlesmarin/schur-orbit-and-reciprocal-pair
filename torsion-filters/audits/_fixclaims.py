# Corrige tres sobre-enunciados del manuscrito y actualiza dos cifras.  Cada uno es un error mio,
# no una cuestion de gusto:
#
# 1. KOSTANT.  El texto decia "the character of an irreducible representation at a REGULAR ELEMENT OF
#    FINITE ORDER lies in {0,+-1}" y lo atribuia a Kostant.  Kostant es sobre el elemento de COXETER
#    / principal, que es lo que dice el propio titulo de Prasad, y que yo mismo verifique esta
#    manana antes de escribirlo mal.  Nuestro elemento NO es el de Coxeter: sus autovalores son
#    xi^{+-1..m} con xi de orden t = h+2.  La Lema no depende de Kostant -- se prueba por el
#    bialternante -- asi que se reescribe como NUESTRA y se cita el mecanismo, no un teorema que no
#    dice eso.
#
# 2. KIM-YACOBI.  Su branching algebra es Sp_{2n} -> Sp_{2n-2}, no Sp_{2R} -> Sp_{2m} x Sp_{2r}, que
#    es la que usamos.  Verificado en su abstract esta misma manana.  Se corrige la atribucion.
#
# 3. LAS CIFRAS.  Se decia "coeficientes 11, 39, 28" leyendo una salida TRUNCADA a 35 caracteres.
#    El volcado completo da termino individual maximo |B tau| = 798 y suma parcial mas lejana 455.
#    La correccion va A FAVOR del enunciado, y por eso mismo habia que hacerla.
p = 'orbit_pair_ii.tex'
s = open(p, encoding='utf-8').read()

CAMBIOS = [
 # 1. Kostant, en la introduccion
 ("""that the character of an irreducible representation at a regular element of finite order lies in
$\\{0,\\pm1\\}$ (Kostant \\cite{Kostant} \\stext; for the current state of that line see \\cite{NPP}), and
that a weight on a wall of the relevant affine Weyl group is""",
  """that a character of a simple module evaluated at a suitable element of finite order is
$0$ or $\\pm1$ --- for the Coxeter element this is Kostant's theorem \\cite{Kostant,Prasad}, and the
line of work it belongs to is active \\cite{NPP}; our element is not the Coxeter element, so we prove
what we need directly. And that a weight on a wall of the relevant affine Weyl group is"""),
 # 1b. la frase de la Lemma que se disculpaba
 ("""Nothing in Lemma~\\ref{lem:T} is new. The dichotomy ``on a wall $\\Rightarrow$ zero, off the walls
$\\Rightarrow$ fold with a sign'' is the standard mechanism of representation theory at roots of
unity, stated for quantum groups in exactly this form by Andersen and Stroppel \\cite{AndersenStroppel}, and the resulting
value lies in $\\{0,\\pm1\\}$ by Kostant's theorem on characters at elements of finite order
\\cite{Kostant}.""",
  """Lemma~\\ref{lem:T} is proved above from the type-$C$ bialternant and owes nothing to what follows,
which is why we can be precise about its standing. Its \\emph{shape} --- on a wall $\\Rightarrow$ zero,
off the walls $\\Rightarrow$ fold with a sign --- is the standard mechanism of representation theory at
roots of unity, stated for quantum groups in exactly this form by Andersen and Stroppel
\\cite{AndersenStroppel}. Its \\emph{species} of conclusion, a character value in $\\{0,\\pm1\\}$ at an
element of finite order, is the species of Kostant's theorem \\cite{Kostant,Prasad}. But that theorem
is about the Coxeter element, and ours is not: the eigenvalues here are $\\xi^{\\pm1},\\dots,\\xi^{\\pm m}$
with $\\xi$ of order $t=h+2$, where $h$ is the Coxeter number. So we neither claim the statement as
new nor inherit it: we prove the case we use."""),
 # 2. Kim-Yacobi
 ("""algebra of Kim and Yacobi \\cite{KimYacobi} \\stext.""",
  """algebra of the rank-one step $\\Sp_{2n}\\downarrow\\Sp_{2n-2}$, for which Kim and Yacobi construct a
standard monomial basis \\cite{KimYacobi} \\stext. Ours is the branching to a \\emph{product},
$\\Sp_{2R}\\downarrow\\Sp_{2m}\\times\\Sp_{2r}$, which is the symmetric-pair setting rather than the
rank-one step; we use the product branching as a computation and do not claim the cited construction
covers it."""),
 # 3. las cifras
 ("""are not small; in a typical case at $t=4$ the surviving $\\eta$ carry coefficients $11$, $39$, $28$
and more, with signs supplied by $\\tau$, and they sum to $\\pm1$.""",
  """are not small. Over the sixteen forms of \\S\\ref{sec:top16} the largest single summand
$|B_{\\eta,\\mu_{\\max}}\\tau_t(\\eta)|$ is $798$, and the running sums travel as far as $455$ from zero
before returning; the filter supplies only signs, and every one of the sixteen lands on $\\pm1$
(Figure~\\ref{fig:collapse})."""),
]

faltan = []
for a, b in CAMBIOS:
    if a in s:
        s = s.replace(a, b, 1)
    else:
        faltan.append(a[:70])

open(p, 'w', encoding='utf-8', newline='\n').write(s)
print('corregidos %d de %d' % (len(CAMBIOS) - len(faltan), len(CAMBIOS)))
for f in faltan:
    print('   NO ENCONTRADO:', repr(f))
