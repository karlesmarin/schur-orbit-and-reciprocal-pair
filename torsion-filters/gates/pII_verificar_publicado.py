# -*- coding: utf-8 -*-
# LO QUE LLEGO CONTRA LO QUE SALIO.  Paper II, arXiv 2608.18302.
#
# La regla de la casa: lo que arXiv publica NO es el fichero que subimos --- es el PDF que arXiv
# compila desde la fuente.  Hay que comparar TEXTO pagina a pagina, no bytes ni tamano.
# Y la referencia no es el PDF del directorio de trabajo (que ya lleva el teorema de la v2): es el
# PDF congelado en el commit del envio.

import hashlib
import subprocess
import sys

import fitz

PUB = r"C:\Users\karles\AppData\Local\Temp\claude\E--proyectos-Curiosity\bbf56496-995b-42b6-9b28-bccf370dd1dc\scratchpad\arxiv_2608.18302.pdf"
ENV = r"C:\Users\karles\AppData\Local\Temp\claude\E--proyectos-Curiosity\bbf56496-995b-42b6-9b28-bccf370dd1dc\scratchpad\enviado_orbit_pair_ii.pdf"
HOY = r"E:\proyectos\Curiosity\research\orbit-pair\paper2\orbit_pair_ii.pdf"


def norm(s):
    return " ".join(s.split())


def paginas(p):
    d = fitz.open(p)
    out = [norm(d[i].get_text()) for i in range(d.page_count)]
    meta = dict(paginas=d.page_count, meta=d.metadata)
    d.close()
    return out, meta


def sha(txt):
    return hashlib.sha256("\n".join(txt).encode("utf-8")).hexdigest()


print("=" * 100)
print("LO QUE LLEGO CONTRA LO QUE SALIO  ---  arXiv 2608.18302 (Paper II)")
print("=" * 100)

pub, mpub = paginas(PUB)
env, menv = paginas(ENV)
hoy, mhoy = paginas(HOY)

print("")
print("  publico arXiv   : %3d paginas   %s" % (mpub["paginas"], PUB.split("\\")[-1]))
print("  enviado (af04581): %3d paginas" % menv["paginas"])
print("  copia de trabajo : %3d paginas   (ya lleva el teorema de la v2)" % mhoy["paginas"])
print("")
print("  SHA-256 del texto extraido")
print("    publico : %s" % sha(pub))
print("    enviado : %s" % sha(env))
print("    trabajo : %s" % sha(hoy))
print("")

if mpub["paginas"] != menv["paginas"]:
    print("  *** EL NUMERO DE PAGINAS NO COINCIDE: %d publico contra %d enviado ***"
          % (mpub["paginas"], menv["paginas"]))

n = min(len(pub), len(env))
malas = []
for i in range(n):
    if pub[i] != env[i]:
        malas.append(i + 1)

print("  COMPARACION PAGINA A PAGINA (publico contra ENVIADO)")
if not malas and len(pub) == len(env):
    print("    %d de %d paginas identicas en texto   -->  LO QUE LLEGO ES LO QUE SALIO" % (n, n))
else:
    print("    paginas que difieren: %d de %d   -->  %s" % (len(malas), n, malas[:40]))
    for i in malas[:3]:
        a, b = pub[i - 1], env[i - 1]
        # primer punto de divergencia
        k = 0
        while k < min(len(a), len(b)) and a[k] == b[k]:
            k += 1
        print("")
        print("    --- pagina %d, primera divergencia en el caracter %d ---" % (i, k))
        print("      publico : ...%s..." % a[max(0, k - 90):k + 120].replace("\n", " "))
        print("      enviado : ...%s..." % b[max(0, k - 90):k + 120].replace("\n", " "))

print("")
print("  Y EL CONTRASTE CON LA COPIA DE TRABAJO (que NO tiene que coincidir: es la v2)")
n2 = min(len(pub), len(hoy))
malas2 = [i + 1 for i in range(n2) if pub[i] != hoy[i]]
print("    paginas que difieren del publico: %d de %d ; paginas %d contra %d"
      % (len(malas2), n2, len(pub), len(hoy)))
print("")
print("  METADATOS del PDF publico")
for k in ("title", "author", "producer", "creator", "creationDate", "modDate"):
    print("    %-13s : %s" % (k, mpub["meta"].get(k)))
print("")
print("=" * 100)
print("DONE")
