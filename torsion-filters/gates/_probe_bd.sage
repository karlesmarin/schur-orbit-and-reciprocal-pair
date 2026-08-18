# sonda: ¿existe la regla de ramificacion B_R' -> B_m' x D_r en esta imagen de Sage?
print("sage ok")
W = WeylCharacterRing("B3")
print("B3 dim (1,0,0) =", W(1,0,0).degree())
try:
    br = branching_rule("B3", "B1xD2", "orthogonal_sum")
    print("regla:", br)
    X = WeylCharacterRing("B1xD2")
    print("branch (1,0,0):", W(1,0,0).branch(X, rule=br))
    print("branch (2,1,0):", W(2,1,0).branch(X, rule=br))
except Exception as e:
    print("FALLA orthogonal_sum:", type(e).__name__, e)
try:
    D = WeylCharacterRing("D2")
    print("D2 (1,0) dim", D(1, 0).degree(), " (1,-1) dim", D(1, -1).degree())
except Exception as e:
    print("FALLA D2:", e)
print("DONE")
