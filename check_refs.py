# -*- coding: utf-8 -*-
"""Two audits a clean LaTeX run cannot do.

(1) TYPE-MATCHED CROSS-REFERENCES.  `\\ref` resolving is not the same as `\\ref` being right: a
    reference reading "Theorem \\ref{lem:foo}" compiles perfectly and is wrong. This maps every label
    to the environment that carries it and flags every citation whose introducing word disagrees.
    (This is the check that caught a table row citing the criterion conjecture when its number came
    from the isolating-witness one.)

(2) BIBLIOGRAPHY IN BOTH DIRECTIONS.  Keys cited but never defined (LaTeX reports these) and keys
    defined but never cited (LaTeX does not).

    python check_refs.py orbit_pair.tex [orbit_pair_es.tex ...]

Authors: Carles Marin, Claude (AI assistant)."""
import io, re, sys

# introducing word -> the environments it may legitimately point at
WORDS = {
    "theorem": {"theorem"}, "teorema": {"theorem"},
    "lemma": {"lemma"}, "lema": {"lemma"},
    "proposition": {"proposition"}, "proposicion": {"proposition"},
    "corollary": {"corollary"}, "corolario": {"corollary"},
    "conjecture": {"conjecture"}, "conjetura": {"conjecture"},
    "question": {"question"}, "pregunta": {"question"},
    "remark": {"remark"}, "observacion": {"remark"},
    "problem": {"problem"}, "problema": {"problem"},
    "example": {"example"}, "ejemplo": {"example"},
    "figure": {"figure"}, "figura": {"figure"},
    "section": {"section", "subsection"}, "seccion": {"section", "subsection"},
}
ABBREV = {"thm": "theorem", "teo": "theorem", "prop": "proposition", "cor": "corollary",
          "conj": "conjecture", "lem": "lemma", "obs": "remark", "fig": "figure"}

NUMBERING = {"theorem", "lemma", "corollary", "proposition", "conjecture", "remark", "question",
             "problem", "example", "definition", "figure", "table", "equation", "align", "gather",
             "multline", "eqnarray"}


def strip_accents(s):
    for a, b in [("ó", "o"), ("ó", "o"), ("é", "e"), ("í", "i"), ("á", "a"), ("ú", "u")]:
        s = s.replace(a, b)
    return s


def audit(path):
    src = io.open(path, encoding="utf-8").read()
    body = src.split("\\begin{thebibliography}")[0]
    print("=" * 78)
    print(path)

    # ---- label -> environment -------------------------------------------------
    # Regexes that try to match \begin{env}[opt]\label{k} in one go break on nested brackets
    # (\begin{theorem}[{\cite[Theorem IX]{LR34}}]) and on braces inside a section title
    # (\section{Proof of Theorem \ref{thm:main}}). Walk the tokens with a stack instead: a label
    # belongs to the innermost open environment, or to the last sectioning command if none is open.
    lab2env = {}
    stack, last_sec = [], None
    tok = re.compile(r"\\begin\{(\w+)\}|\\end\{(\w+)\}|\\(sub)?section\b|\\label\{([^}]+)\}")
    for m in tok.finditer(body):
        beg, end, sub, lab = m.group(1), m.group(2), m.group(3), m.group(4)
        if beg:
            stack.append(beg)
        elif end:
            while stack and stack.pop() != end:
                pass
        elif lab:
            # `document`, `center`, `proof`, `tabular` ... carry no number, so look outward for the
            # innermost environment that actually numbers something; a section label finds none.
            inner = next((e for e in reversed(stack) if e in NUMBERING), None)
            lab2env[lab] = inner or last_sec or "section"
        else:
            last_sec = "subsection" if sub else "section"

    # ---- type-matched references ---------------------------------------------
    bad = 0
    pat = re.compile(r"(\w+)\.?~?\s*\\(?:ref|eqref)\{([^}]+)\}")
    for m in pat.finditer(body):
        word = strip_accents(m.group(1).lower().rstrip("."))
        key = m.group(2)
        if word in ABBREV:
            word = ABBREV[word]
        if word not in WORDS:
            continue
        env = lab2env.get(key)
        if env is None:
            print("  ?? %-28s cited as '%s' but no environment carries that label"
                  % (key, m.group(1)))
            bad += 1
            continue
        if env not in WORDS[word]:
            line = body[:m.start()].count("\n") + 1
            print("  MISMATCH line %-5d '%s \\ref{%s}'  ->  label is a %s"
                  % (line, m.group(1), key, env))
            bad += 1
    print("  type-matched references: %d mismatch(es)" % bad)

    # ---- bibliography both ways ---------------------------------------------
    defined = set(re.findall(r"\\bibitem(?:\[[^\]]*\])?\{([^}]+)\}", src))
    cited = set()
    for m in re.finditer(r"\\cite[a-zA-Z]*(?:\[[^\]]*\])*\{([^}]+)\}", src):
        cited.update(k.strip() for k in m.group(1).split(","))
    miss = sorted(cited - defined)
    orph = sorted(defined - cited)
    print("  bibitems: %d defined, %d cited" % (len(defined), len(cited)))
    if miss:
        print("  CITED BUT NOT DEFINED: %s" % ", ".join(miss))
    if orph:
        print("  DEFINED BUT NEVER CITED: %s" % ", ".join(orph))
    if not miss and not orph:
        print("  bibliography: clean in both directions")
    return bad + len(miss) + len(orph)


total = 0
for p in (sys.argv[1:] or ["orbit_pair.tex"]):
    total += audit(p)
print("=" * 78)
print("TOTAL problems: %d" % total)
