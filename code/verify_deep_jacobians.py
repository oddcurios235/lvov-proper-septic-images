#!/usr/bin/env python3
"""Recompute full original-variable Jacobians at the three q0 rank-one border values."""
import json,pathlib,numpy as np
from verify_deep_minors import ROOT,P,Q0,rank_mod

def border_tuple(B,u):
    N=len(u);n=N+1;A=[]
    for b in B:
        M=np.zeros((n,n),dtype=np.int64);M[:N,:N]=b;A.append(M)
    M=np.zeros((n,n),dtype=np.int64);M[:N,N]=u;A.append(M)
    return A

def evaluate(A):
    n=A[0].shape[0];I=np.eye(n,dtype=np.int64);V=np.zeros((n,n),dtype=np.int64)
    for w,c in Q0.items():
        t=I
        for i in w:t=t@A[i]%P
        V=(V+c*t)%P
    return V

def jacobian(A):
    n=A[0].shape[0];I=np.eye(n,dtype=np.int64)
    J=np.zeros((n*n,7*n*n),dtype=np.int64)
    for w,c in Q0.items():
        L=[I]
        for i in w:L.append(L[-1]@A[i]%P)
        R=[None]*8;R[7]=I
        for j in range(6,-1,-1):R[j]=A[w[j]]@R[j+1]%P
        for j,i in enumerate(w):J[:,i*n*n:(i+1)*n*n]+=c*np.kron(L[j],R[j+1].T)
        J%=P
    return J

def run():
    out=[]
    for c in json.loads((ROOT/'data/deep_core_inputs.json').read_text())['cases']:
        N=c['N'];n=N+1;A=border_tuple([np.array(M,dtype=np.int64)%P for M in c['B']],np.array(c['u'],dtype=np.int64)%P)
        V=evaluate(A);rank=rank_mod(jacobian(A))
        assert rank==n*n-1 and np.any(V[:N,N]) and not np.any(V[:,:N]) and not np.any(V[N,:])
        out.append({'N':N,'matrix_order':n,'full_Jacobian_rank_mod_101':rank,'target_sl_dimension':n*n-1,'rank_one_square_zero_value':bool(not np.any(V@V%P) and np.count_nonzero(V[:N,N])>0)})
        print(out[-1],flush=True)
    (ROOT/'results/deep_jacobians.json').write_text(json.dumps({'prime':P,'cases':out},indent=2)+'\n')

if __name__=='__main__':run()
