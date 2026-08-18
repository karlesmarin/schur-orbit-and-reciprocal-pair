# -*- coding: utf-8 -*-
# Authors: Carles Marin, Claude (AI assistant).
"""CAZA DE CONTRADICCIONES dentro de UNA edicion.

check_parity compara las DOS ediciones; check_numbers ata los numeros a las corridas archivadas.
Ninguno mira si el paper se contradice A SI MISMO: un recuento citado en §8 y otra vez en §9 con
otra cifra, una hipotesis enunciada fuerte en un sitio y debil en otro, una figura cuyo pie dice
algo que el texto niega.

Aqui se agrupan todos los enteros de 2+ cifras del CUERPO por su valor, con su contexto, para que
un par de apariciones del MISMO hecho con cifras distintas salte a la vista.  No decide nada: es un
listado para leer.  Los falsos positivos son la mayoria (un 14 de rango y un 14 de recuento no se
contradicen), y por eso se filtra por sustantivo.
"""
import collections
import io
import re
import sys

NOUNS = ('shape', 'forma', 'case', 'caso', 'fail', 'fallo', 'zero', 'cero', 'witness', 'testigo',
         'collision', 'colision', 'colisi', 'invariant', 'invariante', 'minor', 'menor',
         'configuration', 'configuraci', 'partition', 'partici', 'fibre', 'fibra', 'label',
         'etiqueta', 'check', 'comprobaci', 'disagree', 'discrepan', 'desacuerdo')


def body(path):
    s = io.open(path, encoding='utf-8').read().split('\\begin{thebibliography}')[0]
    if '\\begin{document}' in s:
        s = s.split('\\begin{document}', 1)[1]
    s = re.sub(r'(?<!\\)%.*', '', s)
    return re.sub(r'\s+', ' ', s)


def main(path):
    s = body(path)
    hits = collections.defaultdict(list)
    for m in re.finditer(r'(?<![\d.])(\d{2,})(?![\d])', s):
        v = m.group(1)
        ctx = s[max(0, m.start() - 70):m.end() + 70]
        if any(n in ctx.lower() for n in NOUNS):
            hits[v].append(ctx)
    print('=' * 100)
    print('%s : valores de 2+ cifras que aparecen MAS DE UNA VEZ cerca de un sustantivo de recuento'
          % path)
    print('=' * 100)
    n = 0
    for v in sorted(hits, key=lambda x: (-len(hits[x]), int(x))):
        if len(hits[v]) < 2:
            continue
        n += 1
        print('\n--- %s  (%d apariciones)' % (v, len(hits[v])))
        for c in hits[v]:
            print('    ...%s...' % c.strip())
    print('\ngrupos a revisar: %d' % n)


for p in sys.argv[1:] or ['orbit_pair.tex']:
    main(p)
