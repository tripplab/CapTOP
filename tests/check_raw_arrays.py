#!/usr/bin/env python3
import json, struct, math, sys
from pathlib import Path
cmd=sys.argv[1]; out=Path(sys.argv[2]);
meta=json.loads((out/'captop_grid_metadata.json').read_text())
n=meta['grid']['total_cells']
def doubles(): return list(struct.unpack('<'+'d'*n,(out/'captop_cube_values_f64.raw').read_bytes()))
def u8(): return list((out/'captop_occupied_u8.raw').read_bytes())
def i64(name): return list(struct.unpack('<'+'q'*n,(out/name).read_bytes()))
def eq(a,b):
    assert len(a)==len(b),(a,b)
    for x,y in zip(a,b):
        if isinstance(y,float) and math.isinf(y): assert math.isinf(x) and x>0,(a,b)
        else: assert x==y,(a,b)
if cmd=='single':
    assert (meta['grid']['nx'],meta['grid']['ny'],meta['grid']['nz'])==(1,1,1); assert meta['counts']['occupied_cubes']==1; eq(doubles(),[0.0]); eq(u8(),[1]); eq(i64('captop_element_id_i64.raw'),[10])
elif cmd=='block':
    assert meta['grid']['total_cells']==8 and meta['counts']['occupied_cubes']==8; eq(doubles(),[0.0]*8); eq(u8(),[1]*8)
elif cmd=='sparse':
    assert meta['grid']['total_cells']==8 and meta['counts']['occupied_cubes']==7 and meta['counts']['missing_cubes']==1; assert sum(1 for v in doubles() if math.isinf(v))==1; assert u8().count(0)==1; assert i64('captop_material_i64.raw').count(-1)==1; assert i64('captop_element_id_i64.raw').count(-1)==1
elif cmd=='material': eq(doubles(),[5.0,9.0]); eq(u8(),[1,1]); eq(i64('captop_material_i64.raw'),[5,9])
elif cmd=='material9': eq(doubles(),[math.inf,0.0]); assert meta['counts']['selected_cubes']==1 and meta['counts']['finite_value_cubes']==1
elif cmd=='scalar': eq(doubles(),[1.5,2.5]); eq(u8(),[1,1])
elif cmd=='threshold': eq(doubles(),[math.inf,0.0,0.0]); assert meta['counts']['selected_cubes']==2
