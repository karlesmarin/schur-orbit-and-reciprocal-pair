# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# gor.sage -- independent Gorenstein test: type t = v_p([(A:M)cap O : A]) with M = A cap rad_p.  Compares with 2 sum W = e d.
import sys
def lat(vs):
    return span([vector(QQ,v) for v in vs], ZZ)
def vp(x,p):
    x=abs(Integer(x)); return x.valuation(p) if x else Infinity
def run(m,q,p):
    v=valuation(q,p); qp=q//p**v; E=euler_phi(p**v); D=euler_phi(q)//2; d=euler_phi(qp)//2
    f=1; x=p%qp
    while x not in (1,qp-1):
        x=(x*p)%qp; f+=1
    g=d//f
    K.<z>=CyclotomicField(q)
    a=z+z**-1
    P=a.minpoly()
    Bm=matrix(QQ,[(a**i).list() for i in range(D)])
    def co(y): return Bm.solve_left(vector(QQ,y.list()))
    def mul(y): return matrix(QQ,[co(y*a**i) for i in range(D)])
    R.<T>=K[]
    pol=prod(T-(z**j+z**-j) for j in range(1,m+1))
    Ms=[mul(c) for c in pol.list()[:-1]]
    L=lat([co(K(1))])
    while True:
        rows=list(L.basis())
        L2=lat(rows+[r*M for r in rows for M in Ms])
        if L2==L: break
        L=L2
    A=L; O=lat(identity_matrix(QQ,D).rows())
    detA=A.basis_matrix().det()
    Kp.<w>=CyclotomicField(qp)
    Pq=(w+w**-1).minpoly()
    Fp=PolynomialRing(GF(p),"X")
    assert Fp(P)==Fp(Pq)**E
    psi=Pq(a)
    def radk(i):
        return lat([co(p**s*psi**(i-s)*a**j) for s in range(i+1) for j in range(D)])
    Ainv=[A]
    # conductor
    Cd=O
    for k in range(D):
        Mk=mul(a**k).inverse()
        Cd=Cd.intersection(lat([r*Mk for r in A.basis()]))
    vC=vp(Cd.basis_matrix().det(),p)
    assert vC % d == 0
    e=vC//d
    rads=[radk(i) for i in range(e+2)]
    assert Cd.intersection(lat(rads[e].basis())) == Cd or True
    Ai=[A.intersection(rads[i]) for i in range(e+2)]
    W=[vp(Ai[i+1].basis_matrix().det(),p)-vp(Ai[i].basis_matrix().det(),p) for i in range(e+1)]
    Mx=Ai[1]
    Tt=O
    for b in Mx.basis():
        bb=sum(b[i]*a**i for i in range(D))
        Mb=mul(bb).inverse()
        Tt=Tt.intersection(lat([r*Mb for r in A.basis()]))
    t=vp(detA,p)-vp(Tt.basis_matrix().det(),p)
    lOA=vp(detA,p); lAC=vC-lOA
    pred=(2*sum(W[:e])==e*d); act=(t==1)
    ok = (pred==act) and lOA==sum(d-W[i] for i in range(e)) and lAC==sum(W[:e]) and (e==0 or W[e]==d)
    print("m=%d q=%d p=%d E=%d d=%d g=%d e=%d W=%s lOA=%d lAC=%d type=%d pred_gor=%s gor=%s HS1223=%s %s" % (m,q,p,E,d,g,e,W,lOA,lAC,t,pred,act,(lOA-lAC>=t-1),"OK" if ok else "MISMATCH"))
    sys.stdout.flush()
args=sys.argv[1:]
if args and args[0]=="list":
    QMAX=int(args[1]); DMAX=int(args[2])
    for q in range(5,QMAX+1):
        D=euler_phi(q)//2
        if D>DMAX: continue
        for p in prime_divisors(q):
            v=valuation(q,p); qp=q//p**v
            if qp in (1,2,3,4,6) or euler_phi(p**v)<2: continue
            for m in range(1,(q+1)//2):
                if 2*m>=q or q in (2*m+1,2*m+2): continue
                if (2*m)%qp and (2*m+1)%qp and (2*m+2)%qp: continue
                run(m,q,p)
else:
    run(int(args[0]),int(args[1]),int(args[2]))
