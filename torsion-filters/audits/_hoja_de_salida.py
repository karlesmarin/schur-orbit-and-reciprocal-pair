# -*- coding: utf-8 -*-
# HOJA DE SALIDA DEL PAPER II --- «solo sale verdad confirmada».   19 de agosto de 2026.
#
# POR QUE.  `_claimaudit.py` ya saca cada afirmacion con marcador de estado y los numeros que
# arrastra, «para poder cotejarlos contra los volcados de gates/ uno a uno» --- a mano.  Esto lo
# cierra: para CADA numero del cuerpo del articulo busca que salida archivada lo respalda, y lo que
# no aparece en ninguna sale marcado.  El mismo cepo, pasado antes sobre otro material nuestro,
# saco alli tres cifras sin corrida detras --- dos de ellas correctas, pero sin respaldo, que en
# un articulo enviado es lo mismo que no tenerlas.
#
# QUE NO ES.  No es un verificador semantico: que un numero aparezca en una salida no prueba que sea
# ESE numero. Es un cepo contra el caso real --- la cifra que ya no corresponde a ninguna corrida ---
# y deja una lista corta para mirar a mano.
#
# ALCANCE.  Cuerpo del articulo, de \begin{document} al primer \bibitem, sin comentarios, sin los
# argumentos de \cite/\ref/\label/\includegraphics, y solo numeros de 3+ digitos: por debajo de eso
# la inmensa mayoria son indices, subindices y numeracion, no medidas.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python _hoja_de_salida.py > _hoja_de_salida_OUT.txt 2>&1

import io
import json
import os
import re
import sys
from collections import Counter

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)

TEX = [("inglesa", "orbit_pair_ii.tex"), ("castellana", "orbit_pair_ii_es.tex")]

# Numeros que no son medidas: anos, identificadores de arXiv/DOI, volumenes y paginas de la
# bibliografia que se cuelan por citas en el cuerpo, y los que son parte de un nombre.
IGNORA = set("""1934 1940 1951 1971 1976 1987 1990 1997 1999 2001 2004 2005 2008 2009 2010 2012
2013 2016 2017 2019 2020 2021 2022 2023 2024 2025 2026 105 100 1000""".split())


def normaliza(s):
    # En LaTeX el separador de miles es `\,`, NO la coma: la coma separa elementos de una lista.
    # Meterla aqui pegaba «13, 36, 64, 81, 48» en un solo numero y «August 18, 2026» en otro.
    return re.sub(r"[\s\\]", "", s)


def cuerpo_de(path):
    s = io.open(path, encoding="utf-8", errors="replace").read()
    s = s.split(r"\begin{document}")[-1]
    s = s.split(r"\bibitem")[0]
    s = re.sub(r"(?<!\\)%.*", "", s)                       # comentarios
    s = re.sub(r"\\(cite|ref|eqref|label|pageref|includegraphics)\s*(\[[^\]]*\])?\{[^}]*\}", " ", s)
    s = re.sub(r"\\includegraphics\s*\[[^\]]*\]", " ", s)
    # El separador de miles de LaTeX es `\,` --- barra MAS COMA ---: `190\,443` es UN numero.
    # Quitando solo la barra quedaba la coma y el numero se partia en «190» y «443», dos cifras de
    # tres digitos que el cepo certifica por azar.  Se unen aqui, y SOLO entre digitos, para no
    # tocar las comas que separan elementos de una lista.
    s = re.sub(r"(?<=\d)\\[,;:!\s]+(?=\d)", "", s)
    s = re.sub(r"(?<=\d)~(?=\d)", "", s)
    return s


# ---- indice de evidencia -----------------------------------------------------------------------
EVID = []
for base in [os.path.join(AQUI, "anc"), os.path.join(RAIZ, "gates"), AQUI]:
    if not os.path.isdir(base):
        continue
    for dirpath, _, files in os.walk(base):
        for f in files:
            if not (f.endswith(".txt") or f.endswith(".json") or f.endswith(".md")):
                continue
            # NO indexar la propia salida de esta gate: la corrida anterior escribe en ella las
            # cifras que estaba marcando como SIN FUENTE, y la siguiente las «certifica» contra si
            # misma.  Medido: 780 y 999 pasaron de marcadas a certificadas sin que cambiara nada
            # mas.  Es el artefacto circular de manual, y aqui casi se lo come el instrumento.
            if f.startswith("_hoja_de_salida") or "HOJA_DE_SALIDA" in f:
                continue
            # Ni ninguna nota NUESTRA que comente las cifras: un comentario no es evidencia.  Paso
            # dos veces --- primero con la propia salida, luego con HOJA_DE_SALIDA_HALLAZGOS.md ---,
            # y las dos veces el sintoma fue el mismo: cifras marcadas que pasaban a certificadas sin
            # que ninguna gate hubiera corrido.
            p = os.path.join(dirpath, f)
            try:
                bruto = io.open(p, encoding="utf-8", errors="replace").read()
            except Exception:
                continue
            # Se guardan el texto CRUDO y una version sin separadores.  La busqueda se hace con
            # frontera de digito en las dos: sin eso, «334» casaba dentro de «11334» y el cepo
            # certificaba cualquier cosa --- medido, 99,5% de falsos positivos a tres digitos.
            EVID.append((os.path.relpath(p, RAIZ), bruto, re.sub(r"[\s\\,]", "", bruto)))

# Indice de TOKENS de digitos: cada tirada maximal de digitos que aparece en la evidencia, en el
# texto crudo y en el compactado.  Buscar por token es lo correcto --- «334» no debe casar dentro de
# «11334» --- y ademas convierte la busqueda en O(1) en vez de un regex por fichero.
TOK = re.compile(r"\d+")
INDICE = {}
for f, bruto, plano in EVID:
    for txt in (bruto, plano):
        for m in TOK.finditer(txt):
            INDICE.setdefault(m.group(), set()).add(f)


def certifica(n):
    return sorted(INDICE.get(n, ()))

print("=" * 100)
print("HOJA DE SALIDA DEL PAPER II --- cada cifra del cuerpo, contra su salida archivada")
print("=" * 100)
print("  %d ficheros de evidencia indexados (anc/, ../gates/, paper2/)" % len(EVID))
print("")

NUM = re.compile(r"(?<![\w.])(\d{3,})(?![\w.])")

resumen = {}
for etiqueta, fichero in TEX:
    p = os.path.join(AQUI, fichero)
    if not os.path.exists(p):
        print("  %s : FALTA" % fichero)
        continue
    cuerpo = cuerpo_de(p)
    crudos = [normaliza(m.group()) for m in NUM.finditer(cuerpo)]
    nums = sorted({n for n in crudos if len(n) >= 3 and n not in IGNORA and not n.startswith("0")},
                  key=lambda x: (-len(x), x))
    ok, sin = {}, []
    for n in nums:
        donde = certifica(n)
        if donde:
            ok[n] = donde[:2]
        else:
            sin.append(n)
    print("  == %s : %d cifras distintas de 3+ digitos, %d certificadas, %d SIN FUENTE"
          % (etiqueta, len(nums), len(ok), len(sin)))
    # El desglose por longitud es obligatorio: a tres digitos el cepo tiene ~89% de falso positivo,
    # asi que «certificada» ahi no significa casi nada y no se puede leer junto con las de 5 y 6.
    porlen = Counter(len(n) for n in nums)
    porlen_ok = Counter(len(n) for n in ok)
    print("     por longitud (y lo que vale la certificacion a esa longitud, ver control abajo):")
    for d in sorted(porlen):
        print("        %d digitos : %3d cifras, %3d certificadas" % (d, porlen[d], porlen_ok[d]))
    if sin:
        print("     ---- SIN SALIDA ARCHIVADA ----")
        for n in sin:
            # contexto para poder mirarlo a mano
            i = cuerpo.find(n) if n in cuerpo else -1
            ctx = ""
            for m in NUM.finditer(cuerpo):
                if normaliza(m.group()) == n:
                    a = max(0, m.start() - 90)
                    ctx = " ".join(cuerpo[a:m.end() + 60].split())
                    break
            print("       %-10s %s" % (n, ctx[:150]))
    print("")
    resumen[etiqueta] = dict(total=len(nums), ok=len(ok), sin=sin)

# ---- DERIVADAS: cifras que ninguna gate imprime porque son una cuenta sobre otras que si -------
# No se silencian: se declara la derivacion y se comprueba la aritmetica.  Si la cuenta no cuadra,
# esto falla igual que si la cifra no existiera.
DERIVADAS = [
    ("780", [15, 66, 120, 126, 126, 15, 66, 120, 126], "suma de la columna «weights» de la tabla de (T)"),
    # 2026-08-20.  El cepo lo daba por no certificado y NO lo esta: en gates/regular_in_G_OUT.txt el
    # dual falla t a t --- 3, 24, 102, 340, 875 en los pares hasta t=12 --- y el articulo cita el
    # TOTAL.  La unica aparicion literal de «1344» en gates/ es dentro de «13440», que es otra cosa.
    ("1344", [3, 24, 102, 340, 875],
     "fallos del dual en t PAR, t<=12, de gates/regular_in_G_OUT.txt (en t impar acierta 100%)"),
    # y su compañera de frase, por el mismo motivo: 10+13+72+100+406+575+1908+2724+7823+11178
    ("24809", [10, 13, 72, 100, 406, 575, 1908, 2724, 7823, 11178],
     "pesos probados, 3<=t<=12, misma tabla"),
]
print("  DERIVADAS --- cifras que son una cuenta sobre otras certificadas")
mal_der = []
for valor, sumandos, que in DERIVADAS:
    s = sum(sumandos)
    ok_d = (str(s) == valor)
    print("      %-8s = %s = %d   %s   %s"
          % (valor, " + ".join(str(x) for x in sumandos), s, "OK" if ok_d else "NO CUADRA", que))
    if not ok_d:
        mal_der.append(valor)
print("")

# ---- CONTROL NEGATIVO: .cuanta potencia tiene este cepo? --------------------------------------
# Con 1291 ficheros de evidencia, un numero de TRES cifras aparece en alguno casi por fuerza.  Sin
# medir eso, «115 de 115 certificadas» no dice nada.  Se toman numeros al azar de cada longitud que
# NO estan en el articulo y se mira cuantos «certifica» el cepo: esa es su tasa de falso positivo.
print("  CONTROL NEGATIVO --- tasa de falso positivo del cepo, por numero de digitos")
en_texto = set()
for _, fichero in TEX:
    pp = os.path.join(AQUI, fichero)
    if os.path.exists(pp):
        en_texto |= {normaliza(m.group()) for m in NUM.finditer(cuerpo_de(pp))}
potencia = {}
for d in (3, 4, 5, 6):
    lo, hi = 10 ** (d - 1), 10 ** d
    paso = max(1, (hi - lo) // 250)          # ~250 candidatos por longitud, sin azar
    cand, x = [], lo
    while len(cand) < 200 and x < hi:
        s = str(x)
        if s not in en_texto:
            cand.append(s)
        x += paso
    if not cand:
        continue
    fp = sum(1 for s in cand if certifica(s))
    potencia[d] = (fp, len(cand))
    print("      %d digitos : %3d de %3d numeros AJENOS salen «certificados»  -> falso positivo %5.1f%%   %s"
          % (d, fp, len(cand), 100.0 * fp / len(cand),
             "el cepo NO discrimina a esta longitud" if fp > 0.8 * len(cand) else "discrimina"))
print("")
resumen["control_negativo"] = {str(k): list(v) for k, v in potencia.items()}

# Una derivada DECLARADA Y COMPROBADA arriba deja de estar «sin salida»: su respaldo son los
# sumandos, que si estan archivados.  Sin este descuento el veredicto no podia ponerse verde nunca
# aunque todo estuviera explicado, que es lo que pasaba el 2026-08-20 con 1344 y 24809.  Ojo: solo
# se descuenta lo que CUADRA --- las de `mal_der` siguen contando, y con razon.
_ok_der = {v for v, _, _ in DERIVADAS if v not in mal_der}
for _v in resumen.values():
    if isinstance(_v, dict) and "sin" in _v:
        _v["sin"] = [n for n in _v["sin"] if n not in _ok_der]
tot_sin = sum(len(v["sin"]) for v in resumen.values() if isinstance(v, dict) and "sin" in v)
if _ok_der:
    print("  (descontadas %d derivada(s) declarada(s) y comprobada(s): %s)"
          % (len(_ok_der), ", ".join(sorted(_ok_der))))
    print("")
json.dump(resumen, open(os.path.join(AQUI, "_hoja_de_salida_DUMP.json"), "w"), indent=1)
print("=" * 100)
if tot_sin:
    print("REVISAR: %d cifras sin salida archivada. Correr su gate, o quitarlas del texto." % tot_sin)
else:
    print("TODAS LAS CIFRAS DEL CUERPO TIENEN SALIDA ARCHIVADA DETRAS.")
print("=" * 100)
print("DONE")
