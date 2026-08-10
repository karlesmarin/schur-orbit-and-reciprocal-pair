r"""Does every cross-reference name the right KIND of object, and does every statement get referenced?

`check_refs.py` type-matches references. This is the stricter sibling written for the pre-posting
sweep: it rebuilds, from the source alone, which environment every \label sits in, reads the number
LaTeX gave it from the .aux, and then checks each \ref against the word that introduces it. It also
lists statements that are never referenced, which is not an error but is worth seeing once.

Usage: python check_xref.py orbit_pair.tex [orbit_pair_es.tex]
Authors: Carles Marin, Claude (AI assistant).
"""
import re
import sys
import pathlib

ENVS = ('theorem', 'lemma', 'corollary', 'proposition', 'conjecture',
        'remark', 'example', 'problem', 'question', 'figure')

WORD_EN = {'theorem': 'Theorem', 'lemma': 'Lemma', 'corollary': 'Corollary',
           'proposition': 'Proposition', 'conjecture': 'Conjecture', 'remark': 'Remark',
           'example': 'Example', 'problem': 'Problem', 'question': 'Question',
           'figure': 'Figure'}
WORD_ES = {'theorem': 'Teorema', 'lemma': 'Lema', 'corollary': 'Corolario',
           'proposition': 'Proposición', 'conjecture': 'Conjetura', 'remark': 'Observación',
           'example': 'Ejemplo', 'problem': 'Problema', 'question': 'Pregunta',
           'figure': 'Figura'}
# abbreviations the paper actually uses
ABBREV = {'Thm.': 'theorem', 'Teo.': 'theorem', 'Cor.': 'corollary', 'Lem.': 'lemma',
          'Prop.': 'proposition', 'Conj.': 'conjecture', 'Obs.': 'remark',
          'Rem.': 'remark', 'Fig.': 'figure'}

QUALIFIER = re.compile(
    r'(Theorem|Teorema|Lemma|Lema|Corollary|Corolario|Proposition|Proposici\\?[oó]n|'
    r'Conjecture|Conjetura|Remark|Observaci\\?[oó]n|Example|Ejemplo|Problem|Problema|'
    r'Question|Pregunta|Figure|Figura|Section|Secci\\?[oó]n|'
    r'Thm\.|Teo\.|Cor\.|Lem\.|Prop\.|Conj\.|Obs\.|Rem\.|Fig\.|\\S)?'
    r'[~ ]*\\ref\{([^}]+)\}')


def env_of_labels(tex):
    """Which environment does each \\label sit inside?"""
    out, stack = {}, []
    pat = re.compile(r'\\begin\{(\w+)\}|\\end\{(\w+)\}|\\label\{([^}]+)\}')
    for m in pat.finditer(tex):
        if m.group(1):
            stack.append(m.group(1))
        elif m.group(2):
            if stack and stack[-1] == m.group(2):
                stack.pop()
        elif m.group(3):
            inner = next((e for e in reversed(stack) if e in ENVS), 'section')
            out[m.group(3)] = inner
    return out


def run(path):
    tex = pathlib.Path(path).read_text(encoding='utf-8')
    auxp = pathlib.Path(path).with_suffix('.aux')
    num = dict(re.findall(r'\\newlabel\{([^}]+)\}\{\{([^}]*)\}', auxp.read_text(encoding='utf-8'))) \
        if auxp.exists() else {}
    envof = env_of_labels(tex)
    es = path.endswith('_es.tex')
    word = WORD_ES if es else WORD_EN

    problems, referenced = [], set(re.findall(r'\\eqref\{([^}]+)\}', tex))
    for m in QUALIFIER.finditer(tex):
        q, key = (m.group(1) or '').strip(), m.group(2)
        referenced.add(key)
        if key not in envof:
            problems.append(f'\\ref{{{key}}} points at no \\label in this file')
            continue
        kind = envof[key]
        if not q:
            continue
        if q in ('\\S', 'Section', 'Sección', 'Secci\\\'on'):
            if kind != 'section':
                problems.append(f'"{q} \\ref{{{key}}}" but {key} is a {kind} ({num.get(key, "?")})')
            continue
        want = ABBREV.get(q, None)
        if want is None:
            want = next((e for e, w in word.items() if w.lower() == q.lower().replace('\\', '')), None)
        if want is None:
            continue
        if want != kind:
            problems.append(f'"{q} \\ref{{{key}}}" but {key} is a {kind} ({num.get(key, "?")})')

    statements = {k: v for k, v in envof.items() if v in ENVS and v != 'figure'}
    never = sorted(k for k in statements if k not in referenced)

    print('=' * 78)
    print(path)
    print(f'  labels: {len(envof)}   numbered statements: {len(statements)}   '
          f'references checked: {len(referenced)}')
    print(f'  wrong-kind or dangling references: {len(problems)}')
    for p in problems:
        print('     ' + p)
    print(f'  statements never referenced: {len(never)}'
          + (f' -> {", ".join(f"{k} ({num.get(k, chr(63))})" for k in never)}' if never else ''))
    return len(problems)


if __name__ == '__main__':
    files = sys.argv[1:] or ['orbit_pair.tex', 'orbit_pair_es.tex']
    bad = sum(run(f) for f in files)
    print('=' * 78)
    print(f'TOTAL wrong-kind or dangling references: {bad}')
