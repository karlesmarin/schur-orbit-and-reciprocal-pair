# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# QUE CASILLAS DE UN LOG ESTAN HECHAS.  Es lo que decide que se salta al reanudar un barrido, asi que
# vive en Python y no en grep dentro del .sh (glue-code-decides-so-it-belongs-in-python).
#
# Entiende los DOS formatos de log:
#
#   nuevo (correr_reanudable.sh):  cada casilla termina con una linea  "=== hecha m:q ==="  o
#       "=== fallo m:q rc=N ===".  El bloque se anade al log SOLO al acabar el proceso, asi que un
#       bloque sin marca no existe: la casilla cortada no deja rastro en el log.
#
#   viejo (correr_por_casilla.sh y anteriores):  "--- casilla m:q ---" y, si fallo,
#       "*** casilla m:q TERMINO CON CODIGO N".  Los marcadores pueden estar PEGADOS a mitad de linea
#       (la salida de sage no acaba en salto de linea), por eso se buscan sin anclar.  Un bloque viejo
#       cuenta como hecho SOLO si le sigue otro bloque o "== fin" y no lleva marca de fallo: el ultimo
#       bloque de un log sin "== fin" es una casilla que murio con el proceso, no una hecha.
#
# Uso:
#   python casillas_hechas.py LOG                -> una linea por casilla: "m:q estado"
#   python casillas_hechas.py LOG --saltar       -> sólo las que el runner debe saltar (una por linea)
#   python casillas_hechas.py LOG --saltar --reintentar-fallos   -> las fallidas NO se saltan
#   python casillas_hechas.py LOG --resumen      -> recuento por estado y codigos de fallo
#
# Authors: Carles Marin, Claude (AI assistant).
import re
import sys
from collections import Counter

HECHA = re.compile(r"=== hecha (\S+) ===")
FALLO_NUEVO = re.compile(r"=== fallo (\S+) rc=(\d+) ===")
# \s+ y no " ": con una lista pasada desde PowerShell la celda llevaba "\r" y la marca salia partida
# en dos lineas ("--- casilla 7:132" / " ---").
CAB_VIEJA = re.compile(r"--- casilla (\S+?)\s+---")
FALLO_VIEJO = re.compile(r"\*\*\* casilla (\S+?)\s+TERMINO CON CODIGO (\d+)")
# "== fin: " es el cierre del runner viejo.  NO vale "== fin" a secas: casaria con "== fin de sesion"
# del runner nuevo, y un bloque viejo cortado seguido de una sesion nueva saldria "hecho".
FIN = re.compile(r"== fin: ")
SESION = re.compile(r"== sesion: ")


def estados(texto):
    """dict m:q -> ('hecha', None) | ('fallo', rc) | ('cortada', None).  Lo ultimo visto manda."""
    est = {}
    # formato nuevo
    for mt in HECHA.finditer(texto):
        est[mt.group(1)] = ("hecha", None)
    for mt in FALLO_NUEVO.finditer(texto):
        # un fallo NO pisa una hecha posterior o anterior: si alguna vez termino bien, esta hecha
        if est.get(mt.group(1), ("",))[0] != "hecha":
            est[mt.group(1)] = ("fallo", int(mt.group(2)))
    # formato viejo: trocear por cabeceras
    cabs = list(CAB_VIEJA.finditer(texto))
    for i, mt in enumerate(cabs):
        celda = mt.group(1).strip()
        ini = mt.end()
        fin = cabs[i + 1].start() if i + 1 < len(cabs) else len(texto)
        bloque = texto[ini:fin]
        if HECHA.search(bloque) or FALLO_NUEVO.search(bloque):
            continue  # ya contado en el formato nuevo
        # si entre esta cabecera y la siguiente arranco una sesion nueva, el runner murio DENTRO de
        # esta casilla: la cabecera siguiente no prueba que el proceso de esta terminara.
        ses = SESION.search(bloque)
        if ses:
            bloque = bloque[:ses.start()]
        fv = FALLO_VIEJO.search(bloque)
        cerrado = ((i + 1 < len(cabs)) and not ses) or bool(FIN.search(bloque))
        if fv:
            nuevo = ("fallo", int(fv.group(2)))
        elif cerrado:
            nuevo = ("hecha", None)
        else:
            nuevo = ("cortada", None)
        previo = est.get(celda)
        if previo is None or previo[0] != "hecha":
            est[celda] = nuevo
    return est


def main():
    args = sys.argv[1:]
    if not args:
        print(__doc__ if __doc__ else "uso: casillas_hechas.py LOG [--saltar [--reintentar-fallos] | --resumen]")
        sys.exit(2)
    try:
        with open(args[0], encoding="utf-8", errors="replace") as fh:
            texto = fh.read()
    except FileNotFoundError:
        texto = ""
    est = estados(texto)
    if "--saltar" in args:
        reintentar = "--reintentar-fallos" in args
        # --reintentar-fallos-viejos: solo los fallos del formato viejo.  Esos 137 eran los d^m
        # monomios de L.order (orden_generado.sage); los fallos del runner nuevo ya son con og_orden
        # y reintentarlos en cada relanzamiento quemaria una hora por OOM real.
        viejos = "--reintentar-fallos-viejos" in args
        nuevos_fallidos = set(mt.group(1) for mt in FALLO_NUEVO.finditer(texto))
        for c, (e, _) in est.items():
            if e == "hecha":
                print(c)
            elif e == "fallo":
                es_viejo = c not in nuevos_fallidos
                if not (reintentar or (viejos and es_viejo)):
                    print(c)
    elif "--resumen" in args:
        cnt = Counter(e for (e, _) in est.values())
        rcs = Counter(rc for (e, rc) in est.values() if e == "fallo")
        print("casillas vistas %d   hechas %d   fallidas %d   cortadas %d   codigos %s"
              % (len(est), cnt["hecha"], cnt["fallo"], cnt["cortada"], dict(rcs)))
        for c, (e, rc) in est.items():
            if e != "hecha":
                print("   %-8s %s%s" % (c, e, "" if rc is None else " rc=%d" % rc))
    else:
        for c, (e, rc) in est.items():
            print("%s %s%s" % (c, e, "" if rc is None else " rc=%d" % rc))


if __name__ == "__main__":
    main()
