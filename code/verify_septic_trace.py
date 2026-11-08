#!/usr/bin/env python3
"""Expand proper degree-seven J^3 generators and calculate word, cyclic-trace
and repeated-input restriction ranks for the manuscript."""
import itertools,collections,json,pathlib
from fractions import Fraction
import numpy as np
P=32003
BASE=pathlib.Path(__file__).resolve().parent.parent

def add(a,b,s=1):
 c=a.copy()
 for k,v in b.items():
  c[k]=c.get(k,0)+s*v
  if c[k]==0:del c[k]
 return c

def mul(a,b):
 c=collections.defaultdict(int)
 for w,s in a.items():
  for x,t in b.items():c[w+x]+=s*t
 return {w:v for w,v in c.items() if v}

def atom(i):return {(i,):1}
def bracket(a,b):return add(mul(a,b),mul(b,a),-1)

def matchings(a):
 if not a:yield ();return
 x=a[0]
 for j in a[1:]:
  for v in matchings(tuple(i for i in a if i!=x and i!=j)):yield ((x,j),)+v

def generators():
 for tri in itertools.combinations(range(7),3):
  a,b,c=tri
  lies=[bracket(bracket(atom(a),atom(b)),atom(c)), bracket(bracket(atom(a),atom(c)),atom(b))]
  rest=tuple(i for i in range(7) if i not in tri)
  for li in lies:
   for pairs in matchings(rest):
    x=[li]+[bracket(atom(i),atom(j)) for i,j in pairs]
    for order in itertools.permutations(range(3)):
     y={():1}
     for j in order:y=mul(y,x[j])
     yield y

def cyc(w):return min(w[i:]+w[:i] for i in range(len(w)))

def key(w,shape):
 pos={}
 for i,grp in enumerate(shape):
  for j in grp:pos[j]=i
 return tuple(pos[v] for v in w)

def rank_rows(rows,ncols,prime=P,return_basis=False):
 piv={}; count=0; labels=[]
 for j,row in enumerate(rows):
  v=np.zeros(ncols,dtype=np.int64)
  for k,x in row.items():v[k]=(v[k]+x)%prime
  while True:
   nonzero=np.flatnonzero(v)
   if not len(nonzero):break
   k=int(nonzero[0]);other=piv.get(k)
   if other is None:
    v=(v*pow(int(v[k]),-1,prime))%prime
    piv[k]=v
    labels.append(j);count+=1
    break
   v=(v-int(v[k])*other)%prime
 return ((count,labels,piv) if return_basis else count)

def exact_rational_rank(rows):
 piv={}
 for row in rows:
  v={j:Fraction(c) for j,c in row.items() if c}
  while v:
   k=min(v)
   if k not in piv:
    d=v[k];piv[k]={j:c/d for j,c in v.items()}
    break
   d=v[k]
   for j,c in piv[k].items():
    w=v.get(j,0)-d*c
    if w:v[j]=w
    elif j in v:del v[j]
 return len(piv)

def basis_rows():
 words=list(itertools.permutations(range(7)))
 cs=sorted({cyc(w) for w in words})
 ck={w:i for i,w in enumerate(cs)}
 return ck

def build():
 gs=list(generators());assert len(gs)==1260
 ck=basis_rows();tr=[]
 for g in gs:
  row=collections.defaultdict(int)
  for w,c in g.items():row[ck[cyc(w)]]+=c
  tr.append({k:v for k,v in row.items() if v})
 trace_rank=rank_rows(tr,len(ck))
 exact_trace_rank=exact_rational_rank(tr)
 assert trace_rank==exact_trace_rank==211
 ret={'basis_count':len(gs),'trace_columns':len(ck),'trace_rank_mod_32003':trace_rank,'trace_rank_over_Q':exact_trace_rank,'formal_trace_kernel_dim_if_W_basis':1260-trace_rank}
 print('TRACE',ret,flush=True)
 # Independent full formal-basis verification: raw word rows (5040 columns) sparse pivots.
 wdict={w:i for i,w in enumerate(itertools.permutations(range(7)))}
 def sparse_rank(rows):
  piv={};count=0; peak=0
  for row in rows:
   v={k:(t%P) for k,t in row.items() if t%P}
   while v:
    k=min(v)
    if k not in piv:
     iv=pow(v[k],-1,P);piv[k]={j:(t*iv)%P for j,t in v.items()};count+=1;peak=max(peak,len(v));break
    f=v[k]
    for j,t in piv[k].items():
     q=(v.get(j,0)-f*t)%P
     if q:v[j]=q
     elif j in v:del v[j]
  return count,peak
 wrank,peak=sparse_rank(({wdict[w]:c for w,c in g.items()} for g in gs))
 ret['W_word_rank_mod_32003']=wrank;ret['W_sparse_peak']=peak
 print('W rank',wrank,flush=True)
 assert wrank==1260 and trace_rank==211
 for shapeName,sizes in [('331',[3,3,1]),('2221',[2,2,2,1]),('22111',[2,2,1,1,1]),('211111',[2,1,1,1,1,1])]:
  groups=[];j=0
  for s in sizes:groups.append(list(range(j,j+s)));j+=s
  keys=sorted(set(itertools.permutations(tuple(i for i,s in enumerate(sizes) for _ in range(s)))))
  ix={k:i for i,k in enumerate(keys)}
  rem=[]
  for g in gs:
   row=collections.defaultdict(int)
   for w,c in g.items():row[ix[key(w,groups)]]+=c
   rem.append({k:v for k,v in row.items() if v})
  rankW=rank_rows(rem,len(keys))
  joint=[{**row,**{len(ck)+k:v for k,v in rm.items()}} for row,rm in zip(tr,rem)]
  joint_rank=rank_rows(joint,len(ck)+len(keys))
  ret[shapeName]={'rhoW_rank_mod_32003':rankW,'rhoV_rank_mod_32003':joint_rank-trace_rank,'trace_plus_rho_rank':joint_rank,'rho_outputs':len(keys)}
  print(shapeName,ret[shapeName],flush=True)
 return ret

if __name__=='__main__':
 data=build();o=BASE/'results'/'septic_trace.json';o.write_text(json.dumps(data,indent=2)+'\n');print('SAVED',o)
