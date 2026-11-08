"""Inspect an h5ad file: obs/var columns, obsm/obsp/uns keys, X shape. Read-only, no full load."""
import h5py, sys
import numpy as np

def decode(v):
    if isinstance(v, bytes): return v.decode('utf-8', 'replace')
    return v

path = sys.argv[1]
with h5py.File(path, 'r') as f:
    print('FILE:', path)
    # X shape: stored as CSR sparse (data/indices/indptr) or dense
    if 'X' in f:
        x = f['X']
        if isinstance(x, h5py.Group):
            print('X: sparse group, keys:', list(x.keys()))
            if 'shape' in x.attrs: print('  attrs shape:', x.attrs['shape'])
            if 'indptr' in x: print('  indptr len:', x['indptr'].shape, '-> n_rows =', x['indptr'].shape[0]-1)
            if 'data' in x: print('  nnz:', x['data'].shape[0])
        else:
            print('X: dense, shape:', x.shape, 'dtype:', x.dtype)
    print('layers:', list(f['layers'].keys()) if 'layers' in f else 'none')
    print('obsm:', {k: (v.shape if hasattr(v,'shape') else '?') for k,v in f['obsm'].items()} if 'obsm' in f else 'none')
    print('obsp:', list(f['obsp'].keys()) if 'obsp' in f else 'none')
    print('uns keys:', list(f['uns'].keys()) if 'uns' in f else 'none')
    if 'uns' in f:
        for k in f['uns'].keys():
            try:
                v = f['uns'][k]
                if isinstance(v, h5py.Dataset):
                    print(f'  uns/{k}: dataset shape={v.shape} dtype={v.dtype}')
                else:
                    print(f'  uns/{k}: group keys={list(v.keys())[:10]}')
            except Exception as e:
                print(f'  uns/{k}: <err {e}>')
    obs = f['obs']
    print(f'OBS: n={len(obs["_index"])} rows')
    for k in obs.keys():
        if k == '_index': continue
        ds = obs[k]
        print(f'  obs/{k}: dtype={ds.dtype} shape={ds.shape}', end='')
        try:
            if ds.dtype.kind in 'fiu':
                vals = ds[:1000]
                print(f'  min={np.nanmin(vals):.4g} max={np.nanmax(vals):.4g} mean={np.nanmean(vals):.4g}', end='')
            elif ds.dtype.kind in 'OSU':
                vals = [decode(v) for v in ds[:5]]
                print(f'  examples={vals}', end='')
            elif ds.dtype.kind == 'b':
                print(f'  bool mean={ds[:1000].mean():.3f}', end='')
        except Exception as e:
            print(f'  <read err {e}>', end='')
        # categorical?
        try:
            if 'categories' in ds.attrs or (isinstance(ds, h5py.Dataset) and ds.dtype.kind == 'i' and 'categories' in f['obs'][k].attrs):
                pass
        except Exception: pass
        print()
        # categories subgroup for categorical columns
        try:
            cat = f['obs'][k]
            if isinstance(cat, h5py.Group) and 'categories' in cat:
                cats = [decode(v) for v in cat['categories'][:]]
                print(f'      categories({len(cats)}): {cats[:20]}')
        except Exception: pass
    var = f['var']
    print(f'VAR: n={len(var["_index"])}')
    for k in var.keys():
        if k == '_index': continue
        print(f'  var/{k}: dtype={var[k].dtype}')
