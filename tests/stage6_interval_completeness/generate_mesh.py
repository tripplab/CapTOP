#!/usr/bin/env python3
import sys
case=sys.argv[1]; out=sys.argv[2]
occ=[]
if case=='single': occ=[(0,0,0)]
elif case=='two': occ=[(0,0,0),(2,0,0)]
elif case=='ring': occ=[(i,j,0) for i in range(3) for j in range(3) if not (i==1 and j==1)]
elif case=='shell': occ=[(i,j,k) for i in range(3) for j in range(3) for k in range(3) if not (i==1 and j==1 and k==1)]
else: raise SystemExit(case)
nodes={}; elems=[]
def nid(x,y,z):
    key=(x,y,z)
    if key not in nodes: nodes[key]=len(nodes)+1
    return nodes[key]
for c,(i,j,k) in enumerate(occ,1):
    elems.append((c,[nid(i,j,k),nid(i+1,j,k),nid(i+1,j+1,k),nid(i,j+1,k),nid(i,j,k+1),nid(i+1,j,k+1),nid(i+1,j+1,k+1),nid(i,j+1,k+1)]))
with open(out,'w') as f:
    f.write('MESH "stage6" dimension 3 ElemType Hexahedra Nnode 8\nCoordinates\n')
    inv=sorted(nodes.items(), key=lambda kv: kv[1])
    for (x,y,z),n in inv: f.write(f'{n} {x} {y} {z}\n')
    f.write('End Coordinates\nElements\n')
    for eid,ns in elems: f.write(str(eid)+' '+' '.join(map(str,ns))+' 1\n')
    f.write('End Elements\n')
