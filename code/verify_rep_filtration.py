#!/usr/bin/env python3
"""Compute the S_7 character of the proper-septic trace kernel and
Young-subgroup invariant dimensions from explicit word expansions."""
import itertools,collections,json,math,pathlib,functools
import numpy as np
from verify_septic_trace import generators,cyc,rank_rows,P
BASE=pathlib.Path(__file__).resolve().parent.parent

def partitions(n,lim=None):
 if n==0:yield ();return
 if lim is None:lim=n
 for v in range(min(n,lim),1-1,-1):
  for tail in partitions(n-v,v):yield (v,)+tail

def cycle_type(s):
 seen=set();pieces=[]
 for i in range(len(s)):
  if i not in seen:
   j=i;k=0
   while j not in seen:seen.add(j);k+=1;j=s[j]
   pieces.append(k)
 return tuple(sorted(pieces,reverse=True))

def representative(t):
 a=[];j=0
 for k in t:
  a.extend(list(range(j+1,j+k))+[j]);j+=k
 return tuple(a)

def rim_hook_data(lam,m):
 n=sum(lam)
 if m>n:return []
 out=[]
 for nu in partitions(n-m):
  nn=nu+(0,)*(len(lam)-len(nu))
  if len(nu)>len(lam) or any(a>b for a,b in zip(nn,lam)):continue
  skew={(i,j) for i in range(len(lam)) for j in range(nn[i],lam[i])}
  if not skew and m:return []
  stack=[next(iter(skew))] if skew else [];seen=set(stack)
  while stack:
   i,j=stack.pop()
   for v in ((i+1,j),(i-1,j),(i,j+1),(i,j-1)):
    if v in skew and v not in seen:seen.add(v);stack.append(v)
  if len(seen)!=m:continue
  if any((i+1,j) in skew and (i,j+1) in skew and (i+1,j+1) in skew for i,j in skew):continue
  h=len({i for i,j in skew})-1
  out.append((nu,h))
 return out

@functools.lru_cache(None)
def irreducible_character(lam,mu):
 if sum(lam)!=sum(mu):raise ValueError
 if not mu:return int(not lam)
 return sum((-1)**h*irreducible_character(nu,mu[1:]) for nu,h in rim_hook_data(lam,mu[0]))

def induced_W_char(s):
 # ordered triple, pair, pair; distinguished position of triple (3 copies)
 z=0
 labels=set(range(7))
 for tri in itertools.combinations(range(7),3):
  tr=set(tri)
  if {s[i] for i in tr}!=tr:continue
  rem=labels-tr
  for pa in itertools.combinations(sorted(rem),2):
   pa=set(pa);pb=rem-pa
   if {s[i] for i in pa}!=pa or {s[i] for i in pb}!=pb:continue
   ct=cycle_type(tuple(sorted(tr).index(s[i]) for i in sorted(tr)))
   lie={(1,1,1):2,(2,1):0,(3,):-1}[ct]
   signA=1 if all(s[i]==i for i in pa) else -1
   signB=1 if all(s[i]==i for i in pb) else -1
   z+=lie*signA*signB
 return 3*z

def prod(z):
 a=1
 for x in z:a*=x
 return a

def centralizer_size(mu):
 C=collections.Counter(mu);z=1
 for k,m in C.items():z*=k**m*math.factorial(m)
 return z

def check():
 cs=sorted({cyc(w) for w in itertools.permutations(range(7))});ix={w:j for j,w in enumerate(cs)}
 tr=[]
 for g in generators():
  row=collections.defaultdict(int)
  for w,c in g.items():row[ix[cyc(w)]]+=c
  tr.append({k:v for k,v in row.items() if v})
 rk,_,piv=rank_rows(tr,len(cs),return_basis=True)
 assert rk==211
 keys=sorted(piv); basis=[piv[k] for k in keys]
 types=list(partitions(7)); chars=[]
 for mu in types:
  s=representative(mu); Wchar=induced_W_char(s)
  perm=np.array([ix[cyc(tuple(s[i] for i in w))] for w in cs])
  Tchar=0
  for key,b in zip(keys,basis):
   w=np.zeros(len(cs),dtype=np.int64);w[perm]=b
   diag=0
   for idx in keys:
    if w[idx]:
     t=int(w[idx]);
     if idx==key:diag=t
     w=(w-t*piv[idx])%P
   assert not w.any(),('not invariant',mu)
   Tchar+=diag
  Tchar%=P
  if Tchar>P//2:Tchar-=P
  v=Wchar-Tchar
  chars.append({'cycle_type':mu,'char_W':Wchar,'char_trace_quotient':Tchar,'char_V':v})
  print('CHAR',mu,Wchar,Tchar,v,flush=True)
 dims={};multiplicities={}
 for lam in types:
  dim=irreducible_character(lam,(1,)*7)
  # exact integer by common denom product
  mult=sum(__import__('fractions').Fraction(x['char_V']*irreducible_character(lam,x['cycle_type']),centralizer_size(x['cycle_type'])) for x in chars)
  assert mult.denominator==1,(lam,mult)
  dims[lam]=dim;multiplicities[lam]=int(mult)
 assert sum(dims[k]*multiplicities[k] for k in dims)==1049
 print('multiplicities',[(k,dims[k],multiplicities[k]) for k in types if multiplicities[k]],flush=True)
 filtration={}
 for name,shape in [('331',(3,3,1)),('2221',(2,2,2,1)),('22111',(2,2,1,1,1)),('211111',(2,1,1,1,1,1))]:
  blocks=[];a=0
  for q in shape:blocks.append(tuple(range(a,a+q)));a+=q
  subgroup=[]
  for product in itertools.product(*(itertools.permutations(b) for b in blocks)):
   perm=[None]*7
   for v,b in zip(product,blocks):
    for i,j in zip(b,v):perm[i]=j
   subgroup.append(tuple(perm))
  components=[];K=0
  for lam in types:
   m=multiplicities[lam]
   if not m:continue
   d=sum(irreducible_character(lam,cycle_type(s)) for s in subgroup)
   assert d%len(subgroup)==0,(name,lam,d)
   fix=d//len(subgroup)
   if fix==0:K+=m*dims[lam]
   components.append({'partition':lam,'multiplicity_in_V':m,'degree':dims[lam],'H_fixed_dimension':fix})
  filtration[name]={'common_kernel_dim':K,'H_size':len(subgroup),'components':components}
  print('FILTRATION',name,K,flush=True)
 data={'prime':P,'trace_image_rank':rk,'characters':chars,'V_decomposition':[{'partition':k,'dimension':dims[k],'multiplicity':multiplicities[k]} for k in types if multiplicities[k]],'common_kernel':filtration}
 (BASE/'results'/'representation_filtration.json').write_text(json.dumps(data,indent=2)+'\n')
 assert [filtration[x]['common_kernel_dim'] for x in ('331','2221','22111','211111')]==[762,163,18,0]
 return data

if __name__=='__main__':check()
