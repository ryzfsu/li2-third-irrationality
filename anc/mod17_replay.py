#!/usr/bin/env python3
"""Standard-library replay: rational determinants and native scaled local entries."""
import json
from fractions import Fraction as Q
from math import gcd
from pathlib import Path


def rational_det(matrix):
    A = [[Q(x) for x in row] for row in matrix]
    result = Q(1)
    for k in range(len(A)):
        pivot = next((i for i in range(k,len(A)) if A[i][k]), None)
        if pivot is None:
            return Q(0)
        if pivot != k:
            A[k], A[pivot] = A[pivot], A[k]
            result = -result
        d = A[k][k]
        result *= d
        for i in range(k+1,len(A)):
            z = A[i][k]/d
            for j in range(k+1,len(A)):
                A[i][j] -= z*A[k][j]
    return result


def native_seed(T):
    t, b = [Q(0)], [Q(0)]
    for j in range(1,10):
        t.append(3*t[-1]+Q(1,j))
        b.append(3*b[-1]+Q(1,j*j))
    C, B = [], []
    for i in range(2,6):
        C.append([Q(3**i-3**j,j-i) if i != j else Q(0)
                  for j in range(1,10)])
        B.append([(t[i]-t[j])/Q(j-i) if i != j else b[i]-T*3**i
                  for j in range(1,10)])
    return C+B+[[Q(1,j+1) for j in range(1,10)]],t,b


def p_valuation(x,p):
    assert x
    e=0
    while x % p == 0:
        e+=1
        x//=p
    return e,x


def scaled_quotient(num,den,scale_power,p):
    """p**scale_power * num/den modulo p, where num is already integral."""
    e,unit = p_valuation(den,p)
    assert e <= scale_power
    return 0 if e < scale_power else num*pow(unit % p,-1,p) % p


def reduced_source(p,t):
    rho=p**t
    apow=[1]
    qt,rb=[0],[0]
    for j in range(1,10*rho):
        apow.append(apow[-1]*3 % p)
        qt.append((3*qt[-1]+scaled_quotient(1,j,t,p)) % p)
        rb.append((3*rb[-1]+scaled_quotient(1,j*j,2*t,p)) % p)
    def entry(kind,i,j):
        if kind == "C":
            return (0,0) if i == j else (scaled_quotient(apow[i]-apow[j],j-i,t,p),0)
        if kind == "B":
            return (rb[i],-apow[i] % p) if i == j else (scaled_quotient(qt[i]-qt[j],j-i,t,p),0)
        assert kind == "V"
        return scaled_quotient(1,i+j,t,p),0
    return rho,apow,qt,rb,entry


def qmod(x,p):
    x=Q(x)
    return x.numerator*pow(x.denominator,-1,p) % p


def seed_coeff_entry(kind,i,j,p):
    matrix,t,b=native_seed(Q(0))
    if kind == "C":
        return qmod(matrix[i-2][j-1],p),0
    if kind == "B":
        return qmod(matrix[i+2][j-1],p),(-pow(3,i,p) % p if i==j else 0)
    return pow(j+1,-1,p),0


def modular_det(matrix,p):
    A=[row[:] for row in matrix]
    z=1
    for k in range(len(A)):
        pivot=next((i for i in range(k,len(A)) if A[i][k] % p),None)
        if pivot is None:
            return 0
        if pivot != k:
            A[k],A[pivot]=A[pivot],A[k]
            z=-z
        d=A[k][k] % p
        z=z*d % p
        inv=pow(d,-1,p)
        for i in range(k+1,len(A)):
            if not A[i][k] % p:
                continue
            v=A[i][k]*inv % p
            for j in range(k+1,len(A)):
                A[i][j]=(A[i][j]-v*A[k][j]) % p
            A[i][k]=0
    return z % p


def check_tower(p,t,full=False,det_points=()):
    rho,apow,qt,rb,entry=reduced_source(p,t)
    seed_matrix,seed_t,seed_b=native_seed(Q(0))
    sequence_checks=0
    for j in range(10*rho):
        q,u=divmod(j,rho)
        assert qt[j] == pow(3,u,p)*qmod(seed_t[q],p) % p
        assert rb[j] == pow(3,u,p)*qmod(seed_b[q],p) % p
        sequence_checks+=2
    seed_cache={(kind,i,j):seed_coeff_entry(kind,i,j,p)
                for kind in ("C","B") for i in range(2,6) for j in range(1,10)}
    seed_cache.update({("V",1,j):seed_coeff_entry("V",1,j,p) for j in range(1,10)})
    checks=0
    # Full coefficientwise comparison at the first tower; at the second,
    # all matching blocks plus representatives from every p-adic distance.
    for kind in ("C","B"):
        for i in range(2*rho,6*rho):
            qi,u=divmod(i,rho)
            js=range(rho,10*rho) if full else sorted({
                rho*q+((u+d)%rho) for q in range(1,10)
                for d in [0,1,-1,*[p**e for e in range(t)]]})
            for j in js:
                qj,v=divmod(j,rho)
                expected=tuple(pow(3,u,p)*x % p for x in seed_cache[kind,qi,qj]) if u==v else (0,0)
                assert entry(kind,i,j)==expected, (p,t,kind,i,j,entry(kind,i,j),expected)
                checks+=1
    for ell in range(1,rho+1):
        u=(-ell)%rho
        js=range(rho,10*rho) if full else sorted({
            rho*q+((u+d)%rho) for q in range(1,10)
            for d in [0,1,-1,*[p**e for e in range(t)]]})
        for j in js:
            qj,v=divmod(j,rho)
            expected=seed_cache["V",1,qj] if u==v else (0,0)
            assert entry("V",ell,j)==expected, (p,t,"V",ell,j)
            checks+=1
    determinant_checks=[]
    for W in det_points:
        rows=[(kind,i) for kind in ("C","B") for i in range(2*rho,6*rho)]
        rows += [("V",ell) for ell in range(1,rho+1)]
        A=[]
        for kind,i in rows:
            row=[]
            for j in range(rho,10*rho):
                c0,c1=entry(kind,i,j)
                row.append((c0+W*c1)%p)
            A.append(row)
        actual=modular_det(A,p)
        g=qmod(rational_det(native_seed(Q(W))[0]),p)
        sign=(-1)**(rho*(rho-1)//2)
        expected=sign*pow(3,4*rho*(rho-1),p)*pow(g,rho,p) % p
        assert actual==expected
        determinant_checks.append({"W":W,"native_scaled_determinant_mod_p":actual,"seed_value_mod_p":g,"sign":sign})
    return {"p":p,"t":t,"rho":rho,"coefficient_pair_checks":checks,
            "scaled_sequence_checks":sequence_checks,"full_entry_comparison":full,
            "native_dense_modular_determinant_checks":determinant_checks}


def main():
    seed=json.loads(Path("seed_exact.json").read_text())
    coeffs=list(map(Q,seed["coefficients_ascending"]))
    M,t,b=native_seed(Q(0))
    assert list(map(str,t))==seed["t_j"]
    assert list(map(str,b))==seed["b_j"]
    evaluation_checks=[]
    for T in [Q(0),Q(1),Q(-1),Q(2),Q(1,3),Q(5,7)]:
        value=rational_det(native_seed(T)[0])
        expected=sum(c*T**i for i,c in enumerate(coeffs))
        assert value==expected
        evaluation_checks.append({"T":str(T),"determinant":str(value)})
    towers=[check_tower(11,1,True,[1]),check_tower(17,1,True,[0,1,2]),check_tower(17,2)]
    out={"status":"PASS","seed_independent_fraction_gaussian_checks":evaluation_checks,
         "tower_checks":towers,
         "scope":"Finite tests supplement the symbolic all-t proof; they do not replace it."}
    Path("verification.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))


if __name__ == "__main__":
    main()
