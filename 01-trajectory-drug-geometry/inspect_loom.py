"""Inspect a loom file: matrix shape, row/col attrs, and point-3 feasibility counts."""
import h5py, sys
import numpy as np
from collections import Counter

def dec(v):
    if isinstance(v, bytes): return v.decode('utf-8', 'replace')
    return str(v)

path = sys.argv[1]
with h5py.File(path, 'r') as f:
    print('FILE:', path)
    print('top-level:', list(f.keys()))
    m = f['matrix']
    print('matrix shape:', m.shape, 'dtype:', m.dtype)
    if 'shape' in m.attrs: print('matrix attrs:', dict(m.attrs))
    print('row_attrs:')
    for k in f['row_attrs'].keys():
        ds = f['row_attrs'][k]
        print(f'  {k}: shape={ds.shape} dtype={ds.dtype}', end='')
        try:
            vals = [dec(v) for v in ds[:5]]
            print(f'  examples={vals}', end='')
        except Exception as e: print(f'  <{e}>', end='')
        print()
    print('col_attrs:')
    col = {}
    for k in f['col_attrs'].keys():
        ds = f['col_attrs'][k]
        print(f'  {k}: shape={ds.shape} dtype={ds.dtype}', end='')
        try:
            if ds.dtype.kind in 'fiu':
                vals = ds[:]
                print(f'  min={np.min(vals):.4g} max={np.max(vals):.4g}', end='')
                col[k] = vals
            else:
                vals = np.array([dec(v) for v in ds[:]])
                col[k] = vals
                u, c = np.unique(vals, return_counts=True)
                print(f'  nunique={len(u)}', end='')
                if len(u) <= 25:
                    print(f'  values={dict(zip([dec(v) for v in u], c))}', end='')
                else:
                    print(f'  top5={sorted(zip(c, [dec(v) for v in u]), reverse=True)[:5]}', end='')
        except Exception as e: print(f'  <{e}>', end='')
        print()
    print('layers:', list(f['layers'].keys()) if 'layers' in f else 'none')
    # Point 3 feasibility: non-transformed cells by arm x TissueID
    try:
        ct_key = next((k for k in col if 'cell' in k.lower() or 'type' in k.lower() or 'class' in k.lower()), None)
        tr_key = next((k for k in col if 'treat' in k.lower() or 'drug' in k.lower() or 'condition' in k.lower()), None)
        tid_key = 'TissueID' if 'TissueID' in col else next((k for k in col if 'tissue' in k.lower()), None)
        print(f'--- point3 keys: celltype={ct_key} treatment={tr_key} tissue={tid_key}')
        if ct_key and tr_key and tid_key:
            ct = np.array([dec(v).lower() for v in col[ct_key]])
            tr = np.array([dec(v).lower() for v in col[tr_key]])
            tid = np.array([dec(v) for v in col[tid_key]])
            mal = np.array([('transform' in v or 'malig' in v or 'neoplast' in v or 'tumor' in v or 'glioma' in v) for v in ct])
            print(f'total cells: {len(ct)}  malignant-like: {mal.sum()}  non-transformed: {(~mal).sum()}')
            veh = np.array([('vehicle' in v or 'dmso' in v or 'control' in v) for v in tr])
            print(f'vehicle arm: {veh.sum()}  treated arm: {(~veh).sum()}')
            print('non-transformed by arm:')
            print('  vehicle:', ((~mal) & veh).sum(), ' treated:', ((~mal) & (~veh)).sum())
            print('non-transformed per TissueID x arm:')
            for t in sorted(set(tid)):
                v = ((~mal) & veh & (tid == t)).sum(); d = ((~mal) & (~veh) & (tid == t)).sum()
                print(f'  {t}: vehicle={v} treated={d}')
            # myeloid specifically
            mye = np.array(['myeloid' in v or 'microglia' in v or 'macrophage' in v or 'monocyte' in v for v in ct])
            print(f'myeloid-like cells: {mye.sum()} (vehicle { (mye & veh).sum() }, treated {(mye & (~veh)).sum()})')
    except Exception as e:
        print('point3 check failed:', e)
