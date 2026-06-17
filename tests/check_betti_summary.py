#!/usr/bin/env python3
import json, sys, pathlib
data=json.loads(pathlib.Path(sys.argv[1]).read_text())
for spec in sys.argv[2:]:
    k,v=spec.split('=',1); cur=data
    for part in k.split('.'): cur=cur[part]
    exp=float(v) if '.' in v else int(v)
    if cur!=exp:
        print(f'{k}: expected {exp}, got {cur}', file=sys.stderr); sys.exit(1)
