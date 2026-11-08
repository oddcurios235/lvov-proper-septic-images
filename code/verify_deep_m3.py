#!/usr/bin/env python3
"""Check the frozen M3 matrix-unit evaluation and quadratic trace invariant of q0."""
import json,numpy as np
from verify_deep_minors import ROOT,Q0

def unit(i,j):
    E=np.zeros((3,3),dtype=np.int64);E[i,j]=1;return E

def run():
    positions=[(0,0),(0,1),(0,2),(2,0),(1,0),(0,1),(1,0)]
    A=[unit(*x) for x in positions]
    D=np.zeros((3,3),dtype=np.int64)
    for w,c in Q0.items():
        M=np.eye(3,dtype=np.int64)
        for i in w:M=M@A[i]
        D+=c*M
    value=int(np.trace(D@D))
    assert np.array_equal(D,np.diag([2,-2,0])) and value==8
    # Total cyclic-diagonal grade is zero, so all words on homogeneous inputs remain diagonal.
    grades=[(j-i)%3 for i,j in positions]
    assert sum(grades)%3==0
    out={'matrix_unit_indices_zero_based':positions,'cyclic_diagonal_grades_mod_3':grades,'sum_grades_mod_3':0,'D_q0':D.tolist(),'Q_q0':value}
    (ROOT/'results/deep_m3.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out,flush=True)

if __name__=='__main__':run()
