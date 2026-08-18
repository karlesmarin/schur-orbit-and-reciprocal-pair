# Inserta las dos secciones que faltaban para que el paper sea hermano del companion.
# Se escribe como fichero y no como heredoc: el heredoc del shell se come un nivel de barras, y
# '\\begin' llegaba como '\begin', donde \b es un escape valido (retroceso) y la busqueda fallaba.
p = 'orbit_pair_ii.tex'
s = open(p, encoding='utf-8').read()
verif = open('_sec_verif.tex', encoding='utf-8').read()
openp = open('_sec_open.tex', encoding='utf-8').read()

anchor = r'\section{Attribution}\label{sec:attr}'
assert anchor in s, 'no esta el ancla de Atribucion'
s = s.replace(anchor, verif + '\n' + anchor, 1)

bib = r'\begin{thebibliography}'
assert bib in s, 'no esta la bibliografia'
i = s.index(bib)
s = s[:i] + openp + '\n' + s[i:]

open(p, 'w', encoding='utf-8', newline='\n').write(s)
print('insertadas: Verification (%d lineas) y Open problems (%d lineas)'
      % (verif.count('\n'), openp.count('\n')))
