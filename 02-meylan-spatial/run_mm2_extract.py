"""Runner: execute the MM2 gene extraction R script inside the Colab session."""
import subprocess, sys

r = subprocess.run(
    ["Rscript", "/content/extract_mm2_genes.R"],
    capture_output=True, text=True, timeout=1500,
)
print(r.stdout)
print(r.stderr, file=sys.stderr)
print("Rscript returncode:", r.returncode)
