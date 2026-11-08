# Smoke-test runner for gbm-point1-colab.ipynb (VM, SMOKE=True, 10k Tumor cells).
# Executes the notebook's code cells in order in one Python process, skipping the
# Colab-only shell cells (pip install, gdown download) and running the curl part
# of the CSV cell via subprocess. Prints per-cell timing; stops at first error.
import json, subprocess, sys, time, traceback, os

NB = "/home/hatch/workspace/bioinfo_downloads/gbm-point1/gbm-point1-colab.ipynb"
LOOM_LOCAL = "/home/hatch/workspace/bioinfo_downloads/gbm-point1/b-side/levitin_annotated_gse148842.loom"
OUT_DIR = "/tmp/point1_smoke"
os.makedirs(OUT_DIR, exist_ok=True)

nb = json.load(open(NB))
cells = nb["cells"]
ns = {}
SKIP = {3}   # pip install (already installed in this venv)
FAILED = None

for i, cell in enumerate(cells):
    if cell["cell_type"] != "code":
        continue
    src = "".join(cell["source"]) if isinstance(cell["source"], list) else cell["source"]
    first = src.splitlines()[0] if src.splitlines() else ""
    print(f"\n===== cell {i}: {first[:80]} =====", flush=True)
    t0 = time.time()
    try:
        if i in SKIP:
            print("SKIPPED (Colab-only shell cell)")
            continue
        if i == 1:  # config: exec then override paths for the VM smoke run
            exec(compile(src, f"<cell {i}>", "exec"), ns)
            ns["LOOM_LOCAL"] = LOOM_LOCAL
            ns["SMOKE"] = True
            ns["OUT_DIR"] = OUT_DIR
            os.makedirs(OUT_DIR, exist_ok=True)
            print("overrides: LOOM_LOCAL=local file, SMOKE=True, OUT_DIR=", OUT_DIR)
            continue
        if i == 4:  # gdown download cell: file already present + MD5 verified by md5sum
            print("SKIPPED (loom already on disk, MD5 9ebdd657... verified)")
            continue
        if i == 5:  # CSV download: curl via subprocess, then the python part
            gh = ns["GITHUB_RAW"]
            for name in ["cross_drug_unit_summary.csv", "cross_drug_global_inference.csv"]:
                subprocess.run(["curl", "-sL", "-o", f"{OUT_DIR}/{name}", f"{gh}/{name}"], check=True)
            py = "\n".join(l for l in src.splitlines() if not l.strip().startswith("!"))
            exec(compile(py, f"<cell {i}>", "exec"), ns)
            continue
        exec(compile(src, f"<cell {i}>", "exec"), ns)
    except Exception:
        print(f"!!! CELL {i} FAILED after {time.time()-t0:.1f}s")
        traceback.print_exc()
        FAILED = i
        break
    print(f"--- cell {i} OK ({time.time()-t0:.1f}s) ---", flush=True)

print()
if FAILED is None:
    print("SMOKE TEST PASSED: all code cells executed")
else:
    print(f"SMOKE TEST FAILED at cell {FAILED}")
    sys.exit(1)
