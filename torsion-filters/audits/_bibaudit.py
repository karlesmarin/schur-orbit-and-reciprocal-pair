# Audita la bibliografia del paper II: toda entrada citada, toda cita con entrada, y ningun nombre
# en prosa sin referencia.
#
# NOTA sobre el instrumento: la primera version partia por '\bibitem{' y dejo de ver nada en cuanto
# la bibliografia paso al estilo del companion, '\bibitem[Etiqueta]{clave}'.  Un auditor que no se
# actualiza con el formato dice "0 entradas" y parece que el fichero esta vacio.  Ahora acepta las
# dos formas.
import re

P2 = 'orbit_pair_ii.tex'
s = open(P2, encoding='utf-8').read()
i = s.find(r'\begin{thebibliography}')
j = s.find(r'\end{thebibliography}')
cuerpo, bib = s[:i], s[i:j]

entradas = {}
for m in re.finditer(re.escape('\\bibitem') + r'(?:\[[^\]]*\])?\{([^}]*)\}', bib):
    entradas[m.group(1)] = True

citadas = set()
for m in re.finditer(re.escape('\\cite') + r'\{([^}]*)\}', s):
    for c in m.group(1).split(','):
        citadas.add(c.strip())

print('entradas: %d    claves citadas: %d' % (len(entradas), len(citadas)))
sin_citar = sorted(set(entradas) - citadas)
sin_entrada = sorted(citadas - set(entradas))
print('  entradas nunca citadas : %s' % (sin_citar if sin_citar else 'ninguna'))
print('  citas sin entrada      : %s' % (sin_entrada if sin_entrada else 'ninguna'))

print()
print('nombres en prosa y su cita mas cercana:')
# Los nombres se DERIVAN de la bibliografia, no se escriben a mano.  La version anterior llevaba una
# lista fija y dejo de cubrir, en silencio, todo autor anadido despues -- Pan, Gauss, Zolotarev,
# Lerch, Gross, Landweber, Polo... .  Un auditor con lista fija envejece sin avisar.
APELLIDO = re.compile(r'\b([A-Z][a-z]{2,}(?:\'[A-Z][a-z]+)?)\b')
PALABRAS = {'The', 'Math', 'Amer', 'Monthly', 'Invent', 'Duke', 'Acta', 'Arith', 'Adv', 'Ann',
            'Proc', 'Natl', 'Sci', 'Comm', 'Phys', 'Discrete', 'Finite', 'Fields', 'Appl',
            'Algebraic', 'Combin', 'Graduate', 'Texts', 'Springer', 'Berlin', 'York', 'Dover',
            'Publ', 'Vol', 'Theory', 'Groups', 'Lie', 'Wiley', 'Sons', 'Edition', 'Classical',
            'Introduction', 'Modern', 'Number', 'Preprint', 'Bonner', 'Schriften', 'Invariant',
            'Tableaux', 'Quantum', 'Representations', 'Character', 'Characters', 'Symmetric',
            'Weyl', 'Dynkin', 'Galois', 'Newton', 'Jacobi', 'Schur', 'Dirac', 'Euler',
            'Zolotarev', 'Gauss', 'Frobenius', 'Lerch'}
derivados = set()
# se parte POR bibitem: una ventana de longitud fija con re.S es codiciosa y se come las entradas
# siguientes -- capturaba 16 de 30 y no lo decia.
trozos = bib.split('\\bibitem')[1:]
if len(trozos) != len(entradas):
    raise AssertionError('el troceado ve %d entradas y el conteo %d' % (len(trozos), len(entradas)))
for tr in trozos:
    autores = tr.split('}', 1)[-1].split(r'\emph')[0]
    for a in APELLIDO.findall(autores):
        if a not in PALABRAS:
            derivados.add(a)
NOMBRES = sorted(derivados)
print('  (%d apellidos derivados de la bibliografia)' % len(NOMBRES))
for n in NOMBRES:
    veces = cuerpo.count(n)
    if not veces:
        continue
    # cuantas de esas menciones llevan un \cite a menos de 60 caracteres por detras o por delante
    con = 0
    for m in re.finditer(re.escape(n), cuerpo):
        ventana = cuerpo[max(0, m.start() - 40):m.end() + 60]
        if '\\cite' in ventana:
            con += 1
    marca = '' if con == veces else '   <-- %d sin cita cerca' % (veces - con)
    print('  %-32s %d menciones, %d con cita%s' % (n, veces, con, marca))
