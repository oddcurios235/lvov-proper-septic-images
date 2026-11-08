#!/usr/bin/env python3
"""Verify q0's prescribed integral border-core data modulo 101.
Reports endpoint determinants and the first lexicographically selected full minor
of the normal map for each N; a nonzero residue proves a nonzero integer minor.
"""
import itertools, collections, json, pathlib
import numpy as np

ROOT=pathlib.Path(__file__).resolve().parent.parent
P=101

def atom(i): return {(i,):1}
def add(a,b,sign=1):
    c=collections.defaultdict(int,a)
    for w,v in b.items(): c[w]+=sign*v
    return {w:v for w,v in c.items() if v}
def mul(a,b):
    c=collections.defaultdict(int)
    for w,u in a.items():
        for v,z in b.items(): c[w+v]+=u*z
    return {w:z for w,z in c.items() if z}
def comm(a,b): return add(mul(a,b),mul(b,a),-1)
def q0():
    standard={a:(-1)**sum(a[i]>a[j] for i in range(4) for j in range(i+1,4)) for a in itertools.permutations(range(4))}
    return comm(mul(standard,comm(atom(4),atom(5))),atom(6))
Q0=q0()
assert len(Q0)==96

def matprod(xs,N):
    out=np.eye(N,dtype=np.int64)
    for A in xs:out=out@A%P
    return out

def det_mod(A):
    M=np.array(A,dtype=np.int64,copy=True)%P
    n=M.shape[0];assert M.shape==(n,n);det=1
    for j in range(n):
        opts=np.flatnonzero(M[j:,j])
        if len(opts)==0:return 0
        k=j+int(opts[0]);
        if k!=j:M[[j,k]]=M[[k,j]];det=(-det)%P
        v=int(M[j,j]);det=det*v%P;inv=pow(v,-1,P)
        for k in range(j+1,n):
            if M[k,j]:M[k]=(M[k]-int(M[k,j])*inv*M[j])%P
    return int(det)

def rank_mod(A):
    A=np.array(A,dtype=np.int64,copy=True)%P
    if A.size==0:return 0
    nrow,ncol=A.shape;r=0
    for j in range(ncol):
        if r==nrow:break
        hits=np.flatnonzero(A[r:,j]);
        if not len(hits):continue
        k=r+int(hits[0]);A[[r,k]]=A[[k,r]]
        A[r]=(A[r]*pow(int(A[r,j]),-1,P))%P
        for i in range(r+1,nrow):
            if A[i,j]:A[i]=(A[i]-int(A[i,j])*A[r])%P
        r+=1
    return r

def matrices(B,u):
    N=len(u); assert len(B)==6 and all(x.shape==(N,N) for x in B)
    u=np.asarray(u,dtype=np.int64).reshape(N,1)%P
    R=np.zeros((N,N),dtype=np.int64);Q=np.zeros((N,N),dtype=np.int64)
    terms=[]
    for w,c in Q0.items():
        if w[-1]==6:R=(R+c*matprod([B[k] for k in w[:-1]],N))%P
        if w[0]==6:Q=(Q+c*matprod([B[k] for k in w[1:]],N))%P
        j=w.index(6)
        # In aY X_i b, the T_i trace pairing removes the adjacent X_i.
        # L* instead uses the entire suffix after Y, including X_i.
        near=w[j+1] if j<6 else -1
        terms.append((c,
                      matprod([B[k] for k in w[:j]],N),
                      matprod([B[k] for k in w[j+1:]],N),
                      near,
                      matprod([B[k] for k in w[j+2:]],N) if near>=0 else None))
    cols=[]
    for i in range(N):
        for j in range(N):
            E=np.zeros((N,N),dtype=np.int64);E[i,j]=1
            Ls=np.zeros((N,N),dtype=np.int64)
            T=[np.zeros((N,N),dtype=np.int64) for _ in range(6)]
            for c,a,b,near,tail in terms:
                h=c*b@E@a
                Ls=(Ls+h)%P
                if near>=0:
                    # For T_i, b is the word AFTER the differentiated X_i.
                    T[near]=(T[near]+c*tail@E@a)%P
            cols.append(np.concatenate([Ls.ravel()]+[(x@u).ravel() for x in T])%P)
    return R,Q,np.array(cols,dtype=np.int64).T%P,(R@u).ravel()%P

def selected_rows(G):
    selected=[]
    for i in range(G.shape[0]):
        trial=selected+[i]
        if rank_mod(G[trial,:])>len(selected):selected=trial
        if len(selected)==G.shape[1]:break
    assert len(selected)==G.shape[1]
    return selected

def run():
    cases=json.loads((ROOT/'data/deep_core_inputs.json').read_text())['cases']
    assert [c['N'] for c in cases]==[3,4,5]
    out=[]
    # Independently checked selected-minor residues for the fixed integral inputs.
    expected_minors={3:8,4:91,5:37}
    for c in cases:
        N=c['N'];B=[np.array(M,dtype=np.int64)%P for M in c['B']];u=np.array(c['u'],dtype=np.int64)%P
        R,Q,G,v=matrices(B,u)
        indices=selected_rows(G)
        rdet=det_mod(R);qdet=det_mod(Q);minor=det_mod(G[indices,:]);rank=rank_mod(G)
        assert all([rdet,qdet,minor]) and rank==N*N and np.any(v)
        if minor != expected_minors[N]:
            raise RuntimeError(f'Unexpected normal minor at N={N}: {minor} != {expected_minors[N]}')
        result={'N':N,'matrix_size':N+1,'endpoint_R_det_mod_101':rdet,'endpoint_Q_det_mod_101':qdet,'normal_shape':list(G.shape),'normal_rank_mod_101':rank,'minor_rows_zero_based':indices,'selected_normal_minor_det_mod_101':minor,'border_value_column_mod_101':v.tolist(),'coefficient_degree':N*N+2*N}
        out.append(result);print(result,flush=True)
    assert sum(x['coefficient_degree'] for x in out)==74
    # One mixed-block check: same q0, same free slot, common scalar parameter t.
    left,right=cases[0:2]
    for t in range(1,16):
        B=[]
        for A,C in zip(left['B'],right['B']):
            A=np.array(A,dtype=np.int64);C=np.array(C,dtype=np.int64)
            M=np.zeros((7,7),dtype=np.int64);M[:3,:3]=A;M[3:,3:]=t*C
            B.append(M%P)
        R,Q,G,v=matrices(B,np.array(left['u']+right['u'],dtype=np.int64))
        if det_mod(R) and det_mod(Q) and rank_mod(G)==49:
            transfer={'N':7,'t':t,'R_det_mod_101':det_mod(R),'Q_det_mod_101':det_mod(Q),'normal_rank_mod_101':49}
            break
    else:raise AssertionError('no mixed-block sample among t=1,...,15')
    print('MIXED',transfer,flush=True)
    (ROOT/'results/deep_minors.json').write_text(json.dumps({'prime':P,'cases':out,'sample_transfer':transfer,'degree_Delta_deep':74,'degree_Q':2,'degree_Delta_deep_times_Q':76},indent=2)+'\n')

if __name__=='__main__':run()
