#!/usr/bin/env python3
"""Resume partial GSE330794 files with curl (urllib gets its Range responses
truncated under parallel load; curl serves the same ranges fine).
For each wanted file whose on-disk size != expected, append the missing
tail in 8 MB curl Range requests, verifying each chunk before appending.
"""
import json
import os
import subprocess
import sys

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE330nnn/GSE330794/suppl/GSE330794_RAW.tar"
OUTDIR = os.path.dirname(os.path.abspath(__file__))
CHUNK = 1 * 1024 * 1024  # 1 MB: transfers get cut around 1.3-2 MB, so stay under it


def curl_range(start, end, dest):
    """Fetch tar bytes [start, end] and append to dest. Returns bytes written."""
    cmd = ["curl", "-s", "--retry", "2", "--retry-delay", "3",
           "-r", f"{start}-{end}", BASE]
    out = subprocess.run(cmd, capture_output=True, timeout=600)
    if out.returncode != 0:
        raise RuntimeError(f"curl rc={out.returncode} stderr={out.stderr[:200]!r}")
    if len(out.stdout) != end - start + 1:
        raise RuntimeError(f"short curl chunk: got {len(out.stdout)} want {end - start + 1}")
    with open(dest, "ab") as f:
        f.write(out.stdout)
    return len(out.stdout)


def main():
    idx = {k: tuple(v) for k, v in json.load(open(os.path.join(OUTDIR, "tar_index.json"))).items()}
    want = [l.strip() for l in open(os.path.join(OUTDIR, "download_list.txt")) if l.strip()]
    for name in want:
        dest = os.path.join(OUTDIR, os.path.basename(name))
        data_off, size = idx[name]
        have = os.path.getsize(dest) if os.path.exists(dest) else 0
        if have == size:
            print(f"SKIP {name}", flush=True)
            continue
        print(f"RESUME {name}: have {have} want {size} missing {size - have}", flush=True)
        pos = have
        ok = True
        while pos < size:
            end = min(pos + CHUNK, size) - 1
            done_chunk = False
            for attempt in range(1, 4):
                try:
                    curl_range(data_off + pos, data_off + end, dest)
                    done_chunk = True
                    break
                except Exception as e:
                    print(f"  chunk at {pos} attempt {attempt} failed: {e}", flush=True)
            if not done_chunk:
                ok = False
                break
            pos = os.path.getsize(dest)
        final = os.path.getsize(dest)
        print(f"  {'OK' if final == size else 'FAIL'} {name} final={final} want={size}", flush=True)
        if final != size:
            ok = False
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
