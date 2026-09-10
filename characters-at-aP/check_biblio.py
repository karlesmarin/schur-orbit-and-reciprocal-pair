# -*- coding: utf-8 -*-
# CADA REFERENCIA DE PUBLICACION, CONTRA SU REGISTRO -- Y LAS DOS EDICIONES IGUALES.
# 9 de septiembre de 2026.
#
# POR QUE.  Las otras cinco compuertas miran el CUERPO del articulo: las cifras, las figuras, las
# atribuciones, la estructura.  Ninguna mira la BIBLIOGRAFIA.  El 8-sep esa bibliografia contenia
# ocho afirmaciones sin comprobar (seis articulos citados solo por arXiv y dos con revista pero sin
# volumen ni paginas) y nada en el directorio podia verlas.  Un volumen equivocado resuelve igual
# de bien que uno correcto: el lector llega a OTRO articulo y no se entera.
#
# QUE COMPRUEBA
#   B1  la referencia de revista de cada entrada == la del registro, en las DOS ediciones
#   B2  'revista': null es una afirmacion NEGATIVA: esa entrada no puede llevar cifras de revista
#   B3  PARIDAD: las dos ediciones citan las mismas entradas, con el mismo arXiv y las mismas cifras
#   B4  SENUELOS: un volumen cambiado, una entrada con revista de mas y un arXiv perdido en la
#       castellana TIENEN que saltar.  (El 9-sep la ES no llevaba el arXiv de Kum24 y nadie lo vio.)
#   B5  --online: vuelve a preguntar si algun 'null' ha dejado de serlo -- A LAS DOS FUENTES.
#       La primera version de esta compuerta preguntaba SOLO al journal-ref de arXiv, que lo
#       rellena el autor y falta a menudo: daba por preprints a Kar24 (Proc. Indian Acad. Sci.
#       135, art. 41), KP24 (J. Algebra 676, 69-81) y AF02 (Adv. Appl. Math. 33, 492-511), los
#       tres publicados.  Un 'sin revista' solo vale si lo firman arXiv Y Crossref.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python check_biblio.py > check_biblio_OUT.txt 2>&1
#       python check_biblio.py --online > check_biblio_OUT.txt 2>&1

import io
import json
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

EDICIONES = [("EN", "characters_at_aP.tex"), ("ES", "characters_at_aP_es.tex")]
REG = json.load(io.open("biblio_registro.json", encoding="utf-8"))["entradas"]
fallos = []


def bibitems(texto):
    """{etiqueta: cuerpo} de la bibliografia."""
    i = texto.index(r"\begin{thebibliography}")
    j = texto.index(r"\end{thebibliography}")
    bloque = texto[i:j]
    out = {}
    for m in re.finditer(r"\\bibitem\{([^}]+)\}(.*?)(?=\\bibitem\{|\Z)", bloque, re.S):
        out[m.group(1)] = m.group(2)
    return out


def cola(cuerpo):
    """Lo que queda de una entrada quitando autores, titulo y el arXiv: la revista y nada mas."""
    m = re.search(r"\\emph\{", cuerpo)
    if not m:
        return cuerpo
    k, prof = m.end(), 1
    while k < len(cuerpo) and prof:
        if cuerpo[k] == "{":
            prof += 1
        elif cuerpo[k] == "}":
            prof -= 1
        k += 1
    s = cuerpo[k:]
    s = re.sub(r"\\texttt\{[^}]*\}", " ", s)          # fuera el arXiv
    s = re.sub(r"\([^()]*[Aa]rXiv[^()]*\)", " ", s)   # fuera el aviso de numeracion
    s = s.replace(r"n.\textsuperscript{o}", "no.").replace("~", " ")
    s = re.sub(r"\\[A-Za-z]+", " ", s)
    s = s.replace("{", " ").replace("}", " ")
    return re.sub(r"\s+", " ", s).strip(" ,.;")


def cifras(s):
    return re.findall(r"\d+", s)


def arxiv_de(cuerpo):
    m = re.search(r"arXiv:([\w./-]+?)\}", cuerpo)
    return m.group(1) if m else None


print("=" * 92)
print("B1/B2 -- CADA ENTRADA CONTRA SU REGISTRO, EN LAS DOS EDICIONES")
tex, items = {}, {}
for ed, fn in EDICIONES:
    tex[ed] = io.open(fn, encoding="utf-8").read()
    items[ed] = bibitems(tex[ed])
    print("   %s: %d entradas en la bibliografia" % (ed, len(items[ed])))

for etq, dato in sorted(REG.items()):
    for ed, _ in EDICIONES:
        if etq not in items[ed]:
            print("  FALLA %s/%s: la entrada no existe" % (ed, etq))
            fallos.append("%s/%s ausente" % (ed, etq))
            continue
        cuerpo = items[ed][etq]
        c = cola(cuerpo)
        ax = arxiv_de(cuerpo)
        if ax != dato["arxiv"]:
            print("  FALLA %s/%s: arXiv '%s', el registro dice '%s'" % (ed, etq, ax, dato["arxiv"]))
            fallos.append("%s/%s arXiv" % (ed, etq))
        if dato["revista"] is None:
            if cifras(c):
                print("  FALLA %s/%s: el registro dice SIN REVISTA y la entrada trae cifras: '%s'"
                      % (ed, etq, c))
                fallos.append("%s/%s revista de mas" % (ed, etq))
            else:
                print("   OK   %s/%s  sin revista, y la entrada no la finge" % (ed, etq))
        else:
            esp = re.sub(r"\s+", " ", dato["revista"].replace("~", " ")).strip()
            if cifras(c) != cifras(esp):
                print("  FALLA %s/%s: cifras %s, el registro dice %s   ('%s')"
                      % (ed, etq, cifras(c), cifras(esp), c))
                fallos.append("%s/%s cifras" % (ed, etq))
            else:
                print("   OK   %s/%s  %s" % (ed, etq, esp))

print("=" * 92)
print("B3 -- PARIDAD: las dos ediciones citan lo mismo")
solo_en = sorted(set(items["EN"]) - set(items["ES"]))
solo_es = sorted(set(items["ES"]) - set(items["EN"]))
print(("   OK   " if not solo_en and not solo_es else "  FALLA ")
      + "mismo juego de entradas  (%d vs %d)" % (len(items["EN"]), len(items["ES"]))
      + ("" if not solo_en and not solo_es else "   solo EN=%s  solo ES=%s" % (solo_en, solo_es)))
if solo_en or solo_es:
    fallos.append("juegos de entradas distintos")

desig = []
for etq in sorted(set(items["EN"]) & set(items["ES"])):
    a, b = items["EN"][etq], items["ES"][etq]
    if arxiv_de(a) != arxiv_de(b):
        desig.append("%s: arXiv EN=%s ES=%s" % (etq, arxiv_de(a), arxiv_de(b)))
    elif cifras(cola(a)) != cifras(cola(b)):
        desig.append("%s: cifras EN=%s ES=%s" % (etq, cifras(cola(a)), cifras(cola(b))))
print(("   OK   " if not desig else "  FALLA ") + "mismo arXiv y mismas cifras en cada entrada")
for d in desig:
    print("           " + d)
    fallos.append("paridad " + d.split(":")[0])

print("=" * 92)
print("B4 -- SENUELOS: el control tiene que fallar cuando toca")
def senuelo(nombre, cuerpo_falso, etq, dato, comprueba):
    ok = comprueba(cuerpo_falso, dato)
    print(("   OK   " if ok else "  FALLA ") + "B4: " + nombre + (" SI se detecta" if ok else ""))
    if not ok:
        fallos.append("senuelo " + nombre)

d = REG["CMP06"]
falso = items["EN"]["CMP06"].replace("208 (2007)", "209 (2007)")
senuelo("un volumen cambiado (CMP 208 -> 209)", falso, "CMP06", d,
        lambda c, dd: cifras(cola(c)) != cifras(dd["revista"]))

d = REG["Kar24"]
falso = items["EN"]["Kar24"].replace(r"\texttt{arXiv:2412.17324}.",
                                     r"J.~Algebra 700 (2026), 1--20, \texttt{arXiv:2412.17324}.")
senuelo("una revista inventada donde no hay ninguna (Kar24)", falso, "Kar24", d,
        lambda c, dd: bool(cifras(cola(c))))

d = REG["Kum24"]
falso = re.sub(r"\\texttt\{arXiv:[\w./-]+\}", "", items["ES"]["Kum24"])
senuelo("un arXiv perdido en la castellana (Kum24, el caso real del 9-sep)", falso, "Kum24", d,
        lambda c, dd: arxiv_de(c) != dd["arxiv"])

if "--online" in sys.argv:
    print("=" * 92)
    print("B5 -- ONLINE: preguntar a arXiv Y A CROSSREF si algun 'sin revista' ha dejado de serlo")
    import time
    import urllib.parse
    import urllib.request
    import xml.etree.ElementTree as ET

    def crossref(titulo):
        """La revista de Crossref para ese titulo, o None si no lo tiene."""
        u = ("https://api.crossref.org/works?rows=3&select=title,container-title,volume,page,"
             "published,DOI&query.bibliographic=" + urllib.parse.quote(titulo))
        try:
            items = json.load(urllib.request.urlopen(u, timeout=45))["message"]["items"]
        except Exception:
            return "??"
        pl = lambda s: "".join(c for c in s.lower() if c.isalnum())
        for it in items:
            t = (it.get("title") or [""])[0]
            if pl(t)[:40] and pl(t)[:40] == pl(titulo)[:40]:
                return "%s %s (%s) %s  DOI %s" % (
                    (it.get("container-title") or ["?"])[0], it.get("volume"),
                    (it.get("published", {}).get("date-parts") or [[None]])[0][0],
                    it.get("page") or "", it.get("DOI"))
        return None

    def titulo_de(etq):
        m = re.search(r"\\emph\{(.*?)\}", items["EN"][etq], re.S)
        if not m:
            return ""
        s = re.sub(r"\\[A-Za-z]+", " ", m.group(1))
        s = s.replace("$", " ").replace("{", " ").replace("}", " ").replace("--", "-")
        return re.sub(r"\s+", " ", s).strip()

    ids = [d["arxiv"] for d in REG.values()]
    url = ("http://export.arxiv.org/api/query?id_list=" + ",".join(ids) + "&max_results=100")
    try:
        raw = urllib.request.urlopen(url, timeout=60).read()
        ns = {"a": "http://www.w3.org/2005/Atom", "ax": "http://arxiv.org/schemas/atom"}
        vivo = {}
        for e in ET.fromstring(raw).findall("a:entry", ns):
            aid = re.sub(r"v\d+$", "", e.find("a:id", ns).text.split("/abs/")[1])
            jr = e.find("ax:journal_ref", ns)
            vivo[aid] = re.sub(r"\s+", " ", jr.text).strip() if jr is not None else None
        for etq, dato in sorted(REG.items()):
            hoy = vivo.get(dato["arxiv"], "NO CONTESTA")
            if hoy == "NO CONTESTA":
                print("   ??   %s: arXiv no contesta por %s" % (etq, dato["arxiv"]))
                continue
            if dato["revista"] is not None:
                print("   OK   %s: publicado, el registro lo recoge" % etq)
                continue
            # un 'sin revista' hay que confirmarlo con LAS DOS fuentes
            if hoy is not None:
                print("  FALLA %s: el registro dice SIN REVISTA y arXiv ya dice '%s'" % (etq, hoy))
                fallos.append("registro caducado " + etq)
                continue
            time.sleep(2)
            cr = crossref(titulo_de(etq))
            if cr == "??":
                print("   ??   %s: Crossref no contesta; el 'sin revista' NO queda confirmado" % etq)
            elif cr:
                print("  FALLA %s: arXiv no lo sabe pero Crossref dice '%s'" % (etq, cr))
                fallos.append("registro caducado " + etq)
            else:
                print("   OK   %s: sigue sin revista, firmado por arXiv y por Crossref" % etq)
    except Exception as exc:
        print("   ??   sin red o arXiv caido: %s  (B5 no da veredicto)" % exc)

print("=" * 92)
if fallos:
    print("RESULTADO: %d FALLOS -> %s" % (len(fallos), fallos))
    sys.exit(1)
print("RESULTADO: cada referencia de publicacion esta comprobada contra su fuente,")
print("           las dos ediciones citan lo mismo, y lo que no tiene revista no la finge.")
