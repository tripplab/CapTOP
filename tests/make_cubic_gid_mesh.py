#!/usr/bin/env python3
import argparse, pathlib
ap=argparse.ArgumentParser(); ap.add_argument('output'); ap.add_argument('--nx',type=int,required=True); ap.add_argument('--ny',type=int,required=True); ap.add_argument('--nz',type=int,required=True); ap.add_argument('--occupied',default='all'); ap.add_argument('--spacing',default='1,1,1'); ap.add_argument('--origin',default='0,0,0'); args=ap.parse_args()
sp=[float(x) for x in args.spacing.split(',')]; org=[float(x) for x in args.origin.split(',')]
if args.occupied=='all': occ=[(i,j,k) for k in range(args.nz) for j in range(args.ny) for i in range(args.nx)]
else: occ=[tuple(map(int,c.split(':'))) for c in args.occupied.split(',') if c]
def nid(i,j,k): return 1+i+(args.nx+1)*(j+(args.ny+1)*k)
with open(args.output,'w') as f:
 f.write('Mesh "captop" dimension 3 ElemType Hexahedra Nnode 8\nCoordinates\n')
 for k in range(args.nz+1):
  for j in range(args.ny+1):
   for i in range(args.nx+1):
    f.write(f'{nid(i,j,k)} {org[0]+sp[0]*i:.17g} {org[1]+sp[1]*j:.17g} {org[2]+sp[2]*k:.17g}\n')
 f.write('End Coordinates\nElements\n')
 for e,(i,j,k) in enumerate(occ,1):
  ns=[nid(i,j,k),nid(i+1,j,k),nid(i+1,j+1,k),nid(i,j+1,k),nid(i,j,k+1),nid(i+1,j,k+1),nid(i+1,j+1,k+1),nid(i,j+1,k+1)]
  f.write(str(e)+' '+' '.join(map(str,ns))+' 1\n')
 f.write('End Elements\n')
