# -*- coding: utf-8 -*-
# _control.py -- LIBRERIA, no gate.  Controles que se niegan a certificar cuando no pueden fallar.
#
# POR QUE EXISTE.  El 13 de agosto de 2026, en una sola sesion, escribi DOS controles que no podian
# fallar y no lo vi hasta despues:
#
#   1. closed_form_prof.py, el "señuelo": comparaba H*(beta) contra (F2(beta) == first), o sea las
#      DOS MITADES leidas sobre la variable equivocada.  Cuando las dos se equivocaban a la vez salia
#      "0 desacuerdos" y yo lo lei como que el señuelo pasaba.  Solo al darle un DENOMINADOR --
#      contrastar contra el 'first' MEDIDO -- aparecieron los 10 fallos que tenia que enseñar.
#
#   2. el test de induccion: monte "Phi_r == 0 => Phi_{r-1} == 0" sobre la poblacion de probe(), y
#      probe() devuelve None en cuanto (i) es cierta -- o sea su poblacion es EXACTAMENTE donde
#      Phi != 0.  El antecedente tenia CERO habitantes y la implicacion salia "0 fallos".
#
# Los dos fallos tienen la misma forma y son detectables por la maquina: una tabla de contingencia
# sabe sola si alguna casilla necesaria esta vacia.  Un control cuya casilla contrafactual esta vacia
# NO es un control: es una tautologia con formato de tabla.
#
# COMO SE USA
#
#     from _control import Contingencia
#     c = Contingencia("H* (hipotesis)", "F2 == first (prediccion)")
#     for beta in poblacion:
#         c.add(Hstar(beta), F2(beta) == first(beta), beta)
#     c.informe()                 # imprime tabla, desacuerdos, y AVISOS
#     if not c.valido(): ...      # False si el control no podia fallar
#
# y para implicaciones:
#
#     im = Implicacion("Phi_r == 0", "Phi' == 0")
#     im.add(A, B, beta) ; im.informe()
#
# NO decide si la hipotesis es cierta: decide si la MEDIDA puede decir algo.  Esa es toda su gracia.
#
# Authors: Carles Marin, Claude (AI assistant).
# Self-test: python _control.py   (re-monta los dos fallos de arriba y comprueba que los caza)

from collections import Counter


class Contingencia:
    """Tabla 2x2 entre una hipotesis y una prediccion, que se niega a certificar si no puede fallar."""

    def __init__(self, nombre_h, nombre_p):
        self.h = nombre_h
        self.p = nombre_p
        self.t = Counter()
        self.ej = {}

    def add(self, hip, pred, testigo=None):
        k = (bool(hip), bool(pred))
        self.t[k] += 1
        if testigo is not None and k not in self.ej:
            self.ej[k] = testigo

    # --- lo que hace util a esto -------------------------------------------------------------
    def desacuerdos(self):
        return self.t[(True, False)] + self.t[(False, True)]

    def avisos(self):
        """Lista de razones por las que esta tabla NO puede fallar.  Vacia = el control es portante."""
        a = []
        n = sum(self.t.values())
        if n == 0:
            a.append("la tabla esta VACIA: no se ha medido nada")
            return a
        if self.t[(True, True)] + self.t[(True, False)] == 0:
            a.append("la hipotesis '%s' NO SE CUMPLE NUNCA (%d formas): el acuerdo es vacuo" % (self.h, n))
        if self.t[(False, False)] + self.t[(False, True)] == 0:
            a.append("la hipotesis '%s' se cumple SIEMPRE (%d formas): no discrimina nada" % (self.h, n))
        if self.t[(True, True)] + self.t[(False, True)] == 0:
            a.append("la prediccion '%s' NO se cumple NUNCA: el control no puede confirmar" % self.p)
        if self.t[(False, False)] + self.t[(True, False)] == 0:
            a.append("la prediccion '%s' se cumple SIEMPRE: el control no puede refutar" % self.p)
        if self.desacuerdos() == 0 and self.t[(False, False)] == 0:
            a.append("CONTRAFACTUAL VACIO: no hay ni una forma con hipotesis falsa Y prediccion falsa, "
                     "asi que '0 desacuerdos' no distingue esta hipotesis de una tautologia")
        return a

    def valido(self):
        return not self.avisos()

    def informe(self, indent="     "):
        n = sum(self.t.values())
        print("%s%-34s | %-34s |   n" % (indent, self.h, self.p))
        print("%s%s" % (indent, "-" * 78))
        for k in sorted(self.t, reverse=True):
            print("%s   %-31s | %-34s | %6d%s"
                  % (indent, "SI" if k[0] else "no", "SI" if k[1] else "no", self.t[k],
                     "   ej: %s" % str(self.ej[k])[:34] if k in self.ej else ""))
        d = self.desacuerdos()
        print("%sdesacuerdos: %d de %d" % (indent, d, n))
        for a in self.avisos():
            print("%s*** AVISO: %s" % (indent, a))
        if self.valido():
            print("%sveredicto: %s" % (indent, "EXACTA" if d == 0 else "NO separa (%d desacuerdos)" % d))
        else:
            print("%sveredicto: SIN VALOR -- el control no podia fallar.  Arreglalo antes de leerlo."
                  % indent)
        return self.valido() and d == 0


class Implicacion(Contingencia):
    """A => B, que ademas cuenta los habitantes del ANTECEDENTE antes de dejar leer nada."""

    def informe(self, indent="     "):
        nA = self.t[(True, True)] + self.t[(True, False)]
        print("%shabitantes del antecedente '%s': %d" % (indent, self.h, nA))
        if nA == 0:
            print("%s*** AVISO: antecedente VACIO -- '%s => %s' es vacuamente cierta y no dice nada"
                  % (indent, self.h, self.p))
        fallos = self.t[(True, False)]
        n = sum(self.t.values())
        print("%s%s => %s : %s  (%d fallos sobre %d habitantes, %d formas en total)"
              % (indent, self.h, self.p, "SI" if (fallos == 0 and nA) else "NO", fallos, nA, n))
        for a in self.avisos():
            print("%s*** AVISO: %s" % (indent, a))
        return nA > 0 and fallos == 0


# ===================================================================== AUTOTEST ==================
# Un instrumento que solo se prueba en su mejor caso no vale.  Estos dos casos son MIS DOS FALLOS
# REALES del 13 de agosto, re-montados: la libreria tiene que cazarlos LOS DOS.
if __name__ == "__main__":
    malo = 0
    print("=" * 96)
    print("AUTOTEST -- se re-montan los dos controles que escribi mal, y hay que cazarlos")
    print("=" * 96)
    print("")

    print("  CASO 1  el señuelo de closed_form_prof.py: las dos mitades leidas sobre beta, asi que")
    print("          cuando ambas fallan coinciden y sale '0 desacuerdos'.")
    c = Contingencia("H*(beta) [mal leida]", "F2(beta) == first [mal leida]")
    for _ in range(164):
        c.add(True, True)          # las que iban bien
    for _ in range(8):
        c.add(False, False)        # las 8 en que AMBAS fallan a la vez -> parecen acuerdo
    ok = c.informe()
    caza = not c.valido() or c.t[(False, False)] > 0
    print("          -> desacuerdos = %d, y la casilla no/no tiene %d: el contrafactual EXISTE,"
          % (c.desacuerdos(), c.t[(False, False)]))
    print("             asi que esta tabla SI es portante.  Lo que fallaba en el original era que")
    print("             la 'prediccion' no se contrastaba contra la VERDAD MEDIDA sino contra otra")
    print("             lectura del mismo error.  La libreria no puede detectar eso sola:")
    print("             SE DICE, no se esconde -- hace falta que la prediccion sea independiente.")
    print("")

    print("  CASO 2  el test de induccion vacio: el antecedente no tenia un solo habitante.")
    im = Implicacion("Phi_r == 0", "Phi'(interior) == 0")
    for _ in range(586):
        im.add(False, True)        # A nunca cierto
    leible = im.informe()
    if leible is False:
        print("          -> CAZADO: la libreria se niega a leer la implicacion.")
    else:
        print("          -> *** NO CAZADO ***")
        malo += 1
    print("")

    print("  CASO 3  hipotesis que se cumple SIEMPRE (no discrimina)")
    c3 = Contingencia("hipotesis trivialmente cierta", "prediccion")
    for _ in range(100):
        c3.add(True, True)
    c3.informe()
    if c3.valido():
        print("          -> *** NO CAZADO ***")
        malo += 1
    else:
        print("          -> CAZADO")
    print("")
    print("=" * 96)
    print("AUTOTEST: %s" % ("PASA" if malo == 0 else "FALLA en %d casos" % malo))
    print("")
    print("LIMITE DECLARADO, y es importante: esta libreria caza tablas que no pueden fallar, NO")
    print("predicciones contrastadas contra la fuente equivocada (caso 1).  Para eso no hay atajo")
    print("mecanico: la prediccion tiene que compararse contra una medida INDEPENDIENTE.")
    raise SystemExit(1 if malo else 0)
