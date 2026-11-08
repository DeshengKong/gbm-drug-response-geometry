#!/usr/bin/env python3
"""Build the GSE296577 RDS metadata-inspection notebook for Colab.
Runs pyreadr inside Colab (high-mem) because the local VM OOM-kills on the
3.8 GB RDS (exit 137). Downloads the RDS from NCBI inside Colab (Google
network, no sandbox proxy truncation), inspects meta.data, saves summary.
"""
import nbformat as nbf

NB_PATH = "/home/hatch/workspace/bioinfo_downloads/gse296577/gbm-meylan-rds-inspect.ipynb"
RDS_URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE296nnn/GSE296577/suppl/"
           "GSE296577_GBM_oHSV_seurat_Full_20260107.rds")
RDS_NAME = "GSE296577_GBM_oHSV_seurat_Full_20260107.rds"
SIZE = 3843442210

cells = []

cells.append(nbf.v4.new_markdown_cell(
    "# GSE296577 Meylan Seurat RDS — metadata inspection\n"
    "Local pyreadr OOM-kills on the 3.8 GB RDS (exit 137, 7 GB VM). "
    "Run here on Colab high-mem instead. Downloads the RDS from NCBI "
    "inside Colab, inspects meta.data with pyreadr, saves a summary JSON "
    "and the meta.data table."))

cells.append(nbf.v4.new_code_cell(
    "!pip install -q pyreadr pandas\n"
    "import pyreadr, pandas as pd\n"
    "print('pyreadr', pyreadr.__version__)"))

cells.append(nbf.v4.new_code_cell(
    "import os, subprocess, time\n"
    f"URL = '{RDS_URL}'\n"
    f"DEST = '{RDS_NAME}'\n"
    f"SIZE = {SIZE}\n"
    "t0 = time.time()\n"
    "# resume-capable download; NCBI serves Range, Colab network is direct\n"
    "r = subprocess.run(['curl', '-s', '-L', '-C', '-', '-o', DEST, URL])\n"
    "print('curl rc =', r.returncode)\n"
    "got = os.path.getsize(DEST)\n"
    "print(f'size = {got}  expected = {SIZE}  match = {got == SIZE}  '\n"
    "      f'({time.time()-t0:.0f}s)')\n"
    "assert got == SIZE, 'RDS download incomplete'"))

cells.append(nbf.v4.new_code_cell(
    "import json\n"
    "res = pyreadr.read_r(DEST)\n"
    "print('top-level keys:', list(res.keys()))\n"
    "summary = {'keys': list(res.keys())}\n"
    "for key, obj in res.items():\n"
    "    print(f'\\n=== {key} === {type(obj).__name__}')\n"
    "    if hasattr(obj, 'shape'):\n"
    "        print('  shape:', obj.shape)\n"
    "        summary[key] = {'type': str(type(obj)), 'shape': list(obj.shape)}\n"
    "    if hasattr(obj, 'columns'):\n"
    "        cols = list(obj.columns)\n"
    "        print('  columns:', cols)\n"
    "        summary[key] = {'columns': cols}\n"
    "        for c in cols:\n"
    "            try:\n"
    "                vc = obj[c].value_counts(dropna=False)\n"
    "                if len(vc) <= 40:\n"
    "                    print(f'    {c}: {vc.to_dict()}')\n"
    "                    summary[key][c] = {str(k): int(v) for k, v in vc.items()}\n"
    "                else:\n"
    "                    print(f'    {c}: {len(vc)} distinct values, e.g. {list(vc.index[:5])}')\n"
    "            except Exception as e:\n"
    "                print(f'    {c}: (unreadable: {e})')\n"
    "with open('rds_meta_summary.json', 'w') as f:\n"
    "    json.dump(summary, f, indent=2, default=str)\n"
    "print('\\nsummary written to rds_meta_summary.json')"))

cells.append(nbf.v4.new_code_cell(
    "# Save the meta.data table itself for local use (CSV, gzipped)\n"
    "import glob\n"
    "meta = None\n"
    "for key, obj in res.items():\n"
    "    if hasattr(obj, 'columns') and 'orig.ident' in list(obj.columns):\n"
    "        meta = obj\n"
    "        print('meta.data found under key:', key)\n"
    "        break\n"
    "if meta is None:\n"
    "    # fall back: the largest DataFrame\n"
    "    cands = [(k, v) for k, v in res.items() if hasattr(v, 'shape')]\n"
    "    cands.sort(key=lambda kv: kv[1].shape[0], reverse=True)\n"
    "    meta = cands[0][1]\n"
    "    print('fallback meta under key:', cands[0][0])\n"
    "print('meta shape:', meta.shape)\n"
    "print('survival/response-like columns:',\n"
    "      [c for c in meta.columns if any(s in c.lower() for s in\n"
    "       ['surv', 'pfs', 'os', 'response', 'responder', 'outcome', 'death', 'progress'])])\n"
    "meta.to_csv('meylan_meta.csv.gz', compression='gzip')\n"
    "print('saved meylan_meta.csv.gz')"))

nb = nbf.v4.new_notebook(cells=cells, metadata={
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3.11.0"},
})
nbf.write(nb, NB_PATH)
print("wrote", NB_PATH)
