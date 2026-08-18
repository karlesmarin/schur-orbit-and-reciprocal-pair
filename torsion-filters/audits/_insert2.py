# Inserta las tres secciones nuevas, las filas de verificacion que faltan y las referencias.
# Fichero y no heredoc: el shell se come un nivel de barras.
p = 'orbit_pair_ii.tex'
s = open(p, encoding='utf-8').read()

# --- 1. las tres secciones, antes de "What this settles of the companion paper" -----------------
anc = r'\section{What this settles of the companion paper}\label{sec:answers}'
assert anc in s, 'no esta el ancla de la seccion de respuestas'
nuevas = open('_sec_new.tex', encoding='utf-8').read()
s = s.replace(anc, nuevas + '\n' + anc, 1)

# --- 2. filas nuevas en la tabla de verificacion -------------------------------------------------
anc2 = r'Values in $\{0,\pm1\}$ at a regular torsion element & --- &'
assert anc2 in s, 'no esta el ancla de la tabla'
filas = (
 r'''Prop. \ref{prop:square}: the two factorisations agree & $t=4$, $r=2$ &
$5/5$ coefficient by coefficient; each side rebuilds $\Phi$ & \stverif\\
\quad decoy: the $\varepsilon$ of a different order & same &
disagrees $5/5$ & \stverif\\
Lemma \ref{lem:epsdet}: $\varepsilon$ as a $0/1$ determinant & --- &
proved; and $\varepsilon\in\{0,\pm1\}$ follows & \stproved\\
\eqref{eq:parity}: the parity route against the character route & $t=4$, $r=2$ &
$7/7$ & \stverif\\
\quad decoy: the same without the rearrangement & same &
disagrees $7/7$ & \stverif\\
\eqref{eq:local}: the local form of the $r_i$ & box $\le8$, $R=3$ &
$5368/5368$ under interlacing & \stverif\\
\quad and outside the interlacing & same &
disagrees on $6732$ of $6732$ & \stverif\\
Lemma \ref{lem:checker}: the checkerboard & box $\le8$ &
proved; $3333$ pairs, $0$ failures & \stproved\\
$E^{(t)}$: transfer, triangular, $\{0,\pm1\}$ & $165\times45$ &
density $10.65\%$; $0$ triangularity violations & \stverif\\
\quad M\"obius structure: \emph{absent} & same &
row and column sums span $-4\dots+5$ & \stverif\\
\S\ref{sec:where}: cancellation inside a cell: \emph{no} & $7$ forms &
$5$ of $7$ have singleton cells only & \stverif\\
The atoms at the top weight & $7$ forms &
\textbf{one} or \textbf{three}, signs $2$ vs $1$ & \stverif\\
\quad opposite pair differing by one $t$-ribbon & same &
$2$ of $7$; decoy on other sizes fires in $0$ & \stverif\\
''')
s = s.replace(anc2, filas + anc2, 1)

# --- 3. referencias que faltan -------------------------------------------------------------------
anc3 = r'\bibitem[Yac10]{Yacobi}'
assert anc3 in s, 'no esta el ancla de la bibliografia'
refs = (
 r'''\bibitem[Mol99]{Molev} A.~I.~Molev, \emph{A basis for representations of symplectic Lie
algebras}, Comm. Math. Phys. \textbf{201} (1999), 591--618; and \emph{Gelfand--Tsetlin bases for
classical Lie algebras}, in: Handbook of Algebra 4 (2006); arXiv:math/0211289.

\bibitem[Sun90]{Sundaram} S.~Sundaram, \emph{Tableaux in the representation theory of the classical
Lie groups}, in: Invariant Theory and Tableaux, IMA Vol. Math. Appl. 19, Springer (1990), 191--225.

\bibitem[Wat25]{Watanabe} H.~Watanabe, \emph{Symplectic tableaux and quantum symmetric pairs},
arXiv:2308.01718.

''')
s = s.replace(anc3, refs + anc3, 1)

open(p, 'w', encoding='utf-8', newline='\n').write(s)
print('insertadas: 3 secciones, 13 filas de verificacion, 3 referencias')
