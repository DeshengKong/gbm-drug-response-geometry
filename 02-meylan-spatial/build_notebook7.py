#!/usr/bin/env python3
"""Build the Stage-A structure-probe notebook for the Meylan RDS on Colab.
Runs in the fresh gbm-meylan session (the old gbm-point1 session's kernel is
wedged behind the inspect notebook's nbconvert). Probes how pyreadr
represents the Seurat S4 object: where meta.data, counts, centroids and the
TCR mapping live. Saves a structure report + meta.data CSV if found.
"""
import nbformat as nbf

NB_PATH = "/home/hatch/workspace/bioinfo_downloads/gse296577/gbm-meylan-probe.ipynb"
RDS_NAME = "GSE296577_GBM_oHSV_seurat_Full_20260107.rds"
SIZE = 3843442210

cells = []

cells.append(nbf.v4.new_markdown_cell(
    "# Meylan RDS — Stage A structure probe\n"
    "Fresh session (gbm-meylan). RDS uploaded from the local VM (verified\n"
    "complete there). Goal: learn how pyreadr represents the Seurat S4\n"
    "object — where meta.data, the counts matrix, centroids (images slot)\n"
    "and the TCR clonotype mapping live — before writing the analysis."))

cells.append(nbf.v4.new_code_cell(
    "!pip install -q pyreadr pandas scipy\n"
    "import pyreadr, pandas as pd\n"
    "print('pyreadr', pyreadr.__version__)"))

cells.append(nbf.v4.new_code_cell(
    "import os\n"
    f"DEST = '{RDS_NAME}'\n"
    f"SIZE = {SIZE}\n"
    "got = os.path.getsize(DEST)\n"
    "print(f'size = {got}  expected = {SIZE}  match = {got == SIZE}')\n"
    "assert got == SIZE, 'RDS incomplete'"))

cells.append(nbf.v4.new_code_cell(
    "import json\n"
    "res = pyreadr.read_r(DEST)\n"
    "print('top-level keys:', list(res.keys()))\n"
    "\n"
    "def describe(obj, depth=0, max_depth=3):\n"
    "    t = type(obj).__name__\n"
    "    info = {'type': t}\n"
    "    if hasattr(obj, 'shape'):\n"
    "        info['shape'] = list(obj.shape)\n"
    "    if hasattr(obj, 'columns'):\n"
    "        info['columns'] = list(obj.columns)[:40]\n"
    "    if isinstance(obj, dict):\n"
    "        info['keys'] = list(obj.keys())[:40]\n"
    "        if depth < max_depth:\n"
    "            info['children'] = {k: describe(v, depth + 1, max_depth)\n"
    "                                for k, v in list(obj.items())[:20]}\n"
    "    if hasattr(obj, '__dict__') and not isinstance(obj, (str, bytes)):\n"
    "        try:\n"
    "            d = vars(obj)\n"
    "            if d:\n"
    "                info['attrs'] = list(d.keys())[:40]\n"
    "        except Exception:\n"
    "            pass\n"
    "    return info\n"
    "\n"
    "report = {k: describe(v) for k, v in res.items()}\n"
    "print(json.dumps(report, indent=2, default=str)[:12000])\n"
    "with open('rds_structure.json', 'w') as f:\n"
    "    json.dump(report, f, indent=2, default=str)\n"
    "print('structure report saved')"))

cells.append(nbf.v4.new_code_cell(
    "# Try to locate meta.data (the DataFrame with orig.ident) and save it\n"
    "meta = None\n"
    "meta_key = None\n"
    "for key, obj in res.items():\n"
    "    if hasattr(obj, 'columns') and 'orig.ident' in list(obj.columns):\n"
    "        meta, meta_key = obj, key\n"
    "        break\n"
    "    if isinstance(obj, dict):\n"
    "        for k2, v2 in obj.items():\n"
    "            if hasattr(v2, 'columns') and 'orig.ident' in list(v2.columns):\n"
    "                meta, meta_key = v2, f'{key}.{k2}'\n"
    "                break\n"
    "print('meta.data key:', meta_key)\n"
    "if meta is not None:\n"
    "    print('meta shape:', meta.shape)\n"
    "    print('columns:', list(meta.columns))\n"
    "    meta.to_csv('meylan_meta.csv.gz', compression='gzip')\n"
    "    print('saved meylan_meta.csv.gz')\n"
    "else:\n"
    "    print('meta.data NOT found at top level — see structure report')"))

nb = nbf.v4.new_notebook(cells=cells, metadata={
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3.11.0"},
})
nbf.write(nb, NB_PATH)
print("wrote", NB_PATH)
