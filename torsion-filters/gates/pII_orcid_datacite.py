# -*- coding: utf-8 -*-
# .LLEGA 2608.18302 A ORCID SOLO?  20 de agosto de 2026.
#
# La memoria orcid-is-fed-by-datacite dice: el registro de Carles lo llena DataCite, y una obra
# llega sola SOLO si el DOI lleva su iD en creators[].nameIdentifiers.  arXiv:2608.09619 subio con
# ese campo VACIO porque la cuenta de arXiv no esta vinculada, y hubo que teclearlo a mano.
# Esto comprueba si el Paper II repite el fallo.  Sin suponer nada.

import json
import subprocess
import sys

IDD = "0009-0007-5637-9688"
DOIS = ["10.48550/arXiv.2608.18302", "10.48550/arXiv.2608.09619"]


def get(url, accept="application/json"):
    r = subprocess.run(["curl", "-s", "-H", "Accept: " + accept, "-A", "Mozilla/5.0", url],
                       capture_output=True, text=True, encoding="utf-8", errors="replace")
    return r.stdout


print("=" * 100)
print("ORCID Y DATACITE  ---  .llega el Paper II solo?")
print("=" * 100)
print("")

for doi in DOIS:
    txt = get("https://api.datacite.org/dois/" + doi)
    try:
        d = json.loads(txt)
    except Exception:
        print("  %s  --> DataCite no devuelve JSON (%s)" % (doi, txt[:120].replace("\n", " ")))
        print("")
        continue
    if "data" not in d:
        print("  %s  --> NO ESTA TODAVIA EN DATACITE  (%s)" % (doi, str(d)[:160]))
        print("")
        continue
    a = d["data"]["attributes"]
    print("  %s" % doi)
    print("    estado        : %s   version %s   actualizado %s"
          % (a.get("state"), a.get("version"), a.get("updated")))
    print("    titulo        : %s" % (a.get("titles", [{}])[0].get("title", "")[:90]))
    for c in a.get("creators", []):
        nid = c.get("nameIdentifiers", [])
        print("    creator       : %-28s  nameIdentifiers: %s   %s"
              % (c.get("name", "")[:28], nid if nid else "[] VACIO",
                 "<-- LLEVA EL iD" if nid else "<-- SIN iD: NO LLEGARA SOLO"))
    print("")

# ---- el registro publico de ORCID, con curl y json.load (NUNCA WebFetch: da falsos negativos)
txt = get("https://pub.orcid.org/v3.0/%s/works" % IDD)
try:
    w = json.loads(txt)
except Exception:
    print("  ORCID: no se pudo leer (%s)" % txt[:160])
    sys.exit(0)

grupos = w.get("group", [])
print("  ORCID %s : %d grupos de obras" % (IDD, len(grupos)))
fuentes = {}
enc = []
for g in grupos:
    for s in g.get("work-summary", []):
        fu = (s.get("source", {}).get("source-name", {}) or {}).get("value", "?")
        fuentes[fu] = fuentes.get(fu, 0) + 1
    ids = []
    for e in (g.get("external-ids", {}) or {}).get("external-id", []):
        ids.append(str(e.get("external-id-value", "")))
    if any("2608.18302" in x for x in ids):
        enc.append([(s.get("title", {}).get("title", {}) or {}).get("value", "")[:70]
                    for s in g.get("work-summary", [])])
    if any("2608.09619" in x for x in ids):
        for s in g.get("work-summary", []):
            print("    [2608.09619] put-code %-12s fuente %-22s  %s"
                  % (s.get("put-code"),
                     (s.get("source", {}).get("source-name", {}) or {}).get("value", "?"),
                     (s.get("title", {}).get("title", {}) or {}).get("value", "")[:60]))
print("    fuentes : %s" % fuentes)
print("")
print("    .esta 2608.18302 (Paper II) en ORCID? : %s" % ("SI -- %s" % enc if enc else "NO"))
print("")
print("=" * 100)
print("DONE")
