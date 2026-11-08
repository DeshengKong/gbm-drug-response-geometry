#!/usr/bin/env python3
"""Sequential curl resume for the GSE296577 RDS file.
The parallel urllib downloader gets its Range responses truncated once the
proxy throttles sustained parallel load (hundreds of round-1 chunks short by
a few KB, remainder fetches returning nothing). Sequential 1 MB curl ranges
proved reliable for the GSE330794 resume.
Scans the file for all-zero / missing 1 MB chunks and fills them in order.
Usage: resume_rds_curl.py
"""
import os
import subprocess
import time

URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE296nnn/GSE296577/suppl/GSE296577_GBM_oHSV_seurat_Full_20260107.rds"
DEST = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "GSE296577_GBM_oHSV_seurat_Full_20260107.rds")
SIZE = 3843442210
CHUNK = 1 * 1024 * 1024


def chunk_is_zero(fd, start, length):
    os.lseek(fd, start, os.SEEK_SET)
    data = os.read(fd, length)
    return len(data) < length or data == b"\x00" * len(data)


def fetch_chunk(start, end):
    for attempt in range(1, 6):
        cmd = ["curl", "-s", "-r", f"{start}-{end}", URL]
        out = subprocess.run(cmd, capture_output=True, timeout=300)
        if out.returncode == 0 and len(out.stdout) == end - start + 1:
            return out.stdout
        time.sleep(2 * attempt)
    raise RuntimeError(f"chunk {start} failed after 5 attempts")


def main():
    n_chunks = (SIZE + CHUNK - 1) // CHUNK
    fd = os.open(DEST, os.O_RDWR)
    missing = []
    for i in range(n_chunks):
        start = i * CHUNK
        end = min(start + CHUNK, SIZE) - 1
        if chunk_is_zero(fd, start, end - start + 1):
            missing.append((start, end))
    print(f"missing chunks: {len(missing)}/{n_chunks}", flush=True)
    for j, (start, end) in enumerate(missing):
        data = fetch_chunk(start, end)
        os.pwrite(fd, data, start)
        if (j + 1) % 50 == 0:
            print(f"  filled {j + 1}/{len(missing)}", flush=True)
        time.sleep(0.2)
    os.close(fd)
    final = os.path.getsize(DEST)
    print(f"DONE size={final} expected={SIZE} match={final == SIZE}", flush=True)


if __name__ == "__main__":
    main()
