#!/usr/bin/env python3
"""Extract selected files from the GSE330794 RAW.tar via HTTP Range requests.
The individual sample files are NOT served as separate HTTP objects by NCBI GEO;
they only exist inside the 213 GB tar. This script scans tar headers (512 B each)
to build an index, then downloads only the wanted byte ranges (~3.6 GB).
"""
import json
import os
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE330nnn/GSE330794/suppl/GSE330794_RAW.tar"
OUTDIR = os.path.dirname(os.path.abspath(__file__))
INDEX_JSON = os.path.join(OUTDIR, "tar_index.json")
WANT_TXT = os.path.join(OUTDIR, "download_list.txt")
LOG = os.path.join(OUTDIR, "extract.log")


def fetch_range(start, end):
    req = urllib.request.Request(BASE, headers={"Range": f"bytes={start}-{end}"})
    with urllib.request.urlopen(req, timeout=180) as r:
        if r.status != 206:
            raise RuntimeError(f"range {start}-{end} -> HTTP {r.status}")
        return r.read()


def scan():
    index = {}
    off = 0
    pending_long = None
    while True:
        hdr = fetch_range(off, off + 511)
        if hdr == b"\x00" * 512:
            break
        name = hdr[0:100].split(b"\x00")[0].decode("utf-8", "replace").lstrip("/")
        prefix = hdr[345:500].split(b"\x00")[0].decode("utf-8", "replace")
        if prefix:
            name = prefix.rstrip("/") + "/" + name
        size = int(hdr[124:136].split(b"\x00")[0].strip() or b"0", 8)
        typeflag = chr(hdr[156]) if hdr[156] else "0"
        data_off = off + 512
        if typeflag == "L":
            raw = fetch_range(data_off, data_off + size - 1) if size else b""
            pending_long = raw.split(b"\x00")[0].decode("utf-8", "replace").lstrip("/")
        elif typeflag in ("x", "g"):
            pass  # pax extended header: skip
        elif typeflag in ("0", "\x00"):
            if pending_long:
                name = pending_long
                pending_long = None
            index[name] = (data_off, size)
        else:
            print(f"WARN unknown typeflag {typeflag!r} at {off} name {name}", flush=True)
        off = data_off + ((size + 511) // 512) * 512
    return index


import http.client

CHUNK = 1 * 1024 * 1024  # 1 MB ranges: less likely to be truncated
PIECE = 256 * 1024


def fetch_range_robust(start, end):
    """Return exactly end-start+1 bytes, resuming the remainder on truncation."""
    buf = bytearray()
    pos = start
    for _round in range(6):
        if pos > end:
            break
        req = urllib.request.Request(BASE, headers={"Range": f"bytes={pos}-{end}"})
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                if r.status != 206:
                    raise RuntimeError(f"HTTP {r.status}")
                while pos <= end:
                    try:
                        piece = r.read(min(PIECE, end - pos + 1))
                    except http.client.IncompleteRead as e:
                        piece = e.partial
                    if not piece:
                        break
                    buf += piece
                    pos += len(piece)
        except Exception:
            pass  # retry the remainder in the next round
        if pos > end:
            break
    if len(buf) != end - start + 1:
        raise RuntimeError(f"short after retries: got {len(buf)} want {end - start + 1}")
    return bytes(buf)


def download_one(args):
    name, (data_off, size) = args
    dest = os.path.join(OUTDIR, os.path.basename(name))
    if os.path.exists(dest) and os.path.getsize(dest) == size:
        return f"SKIP {name}"
    for attempt in range(1, 4):
        try:
            # resume from the last whole-chunk boundary; drop a possibly partial tail chunk
            done = os.path.getsize(dest) if os.path.exists(dest) else 0
            done = (done // CHUNK) * CHUNK
            if os.path.exists(dest) and os.path.getsize(dest) != done:
                with open(dest, "r+b") as f:
                    f.truncate(done)
            with open(dest, "ab" if done else "wb") as f:
                pos = done
                while pos < size:
                    end = min(pos + CHUNK, size) - 1
                    data = fetch_range_robust(data_off + pos, data_off + end)
                    f.write(data)
                    f.flush()
                    pos += len(data)
            if os.path.getsize(dest) == size:
                return f"OK {name}"
            raise RuntimeError("size mismatch after chunks")
        except Exception as e:
            if attempt == 3:
                return f"FAIL {name}: {e}"
    return f"FAIL {name}"


def main():
    if os.path.exists(INDEX_JSON):
        index = {k: tuple(v) for k, v in json.load(open(INDEX_JSON)).items()}
        print(f"loaded index: {len(index)} entries", flush=True)
    else:
        print("scanning tar headers ...", flush=True)
        index = scan()
        json.dump(index, open(INDEX_JSON, "w"))
        print(f"scanned index: {len(index)} entries", flush=True)
    want = [l.strip() for l in open(WANT_TXT) if l.strip()]
    missing = [w for w in want if w not in index]
    if missing:
        print("MISSING from tar index:", missing, flush=True)
        sys.exit(1)
    jobs = [(w, index[w]) for w in want]
    total = sum(index[w][1] for w in want)
    print(f"downloading {len(jobs)} files, {total/1e9:.2f} GB", flush=True)
    with ThreadPoolExecutor(max_workers=6) as ex:
        for res in ex.map(download_one, jobs):
            print(res, flush=True)
            with open(LOG, "a") as lf:
                lf.write(res + "\n")
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
