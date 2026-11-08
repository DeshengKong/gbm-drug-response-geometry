#!/usr/bin/env python3
"""Robust parallel chunked Range download for a single large file.
- Each chunk is fetched with resume-within-chunk: if a Range response is
  truncated mid-stream, the remainder is re-requested and appended.
- A chunk that still fails after retries is reported, not fatal: the main
  loop re-queues failed chunks for up to 3 rounds.
- Resumes from the largest whole-chunk boundary of an existing sequential prefix.
Usage: parallel_range_download.py <url> <dest> <size_bytes>
"""
import http.client
import os
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor

URL, DEST, SIZE = sys.argv[1], sys.argv[2], int(sys.argv[3])
CHUNK = 1 * 1024 * 1024  # 1 MB ranges: less likely to be truncated
WORKERS = 6
PIECE = 256 * 1024


def fetch_range_robust(start, end):
    """Return exactly end-start+1 bytes, resuming the remainder on truncation."""
    buf = bytearray()
    pos = start
    for _round in range(6):
        if pos > end:
            break
        req = urllib.request.Request(URL, headers={"Range": f"bytes={pos}-{end}"})
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


def job(args):
    start, end = args
    try:
        data = fetch_range_robust(start, end)
        fd = os.open(DEST, os.O_WRONLY)
        try:
            os.pwrite(fd, data, start)
        finally:
            os.close(fd)
        return (start, True, "")
    except Exception as e:
        return (start, False, str(e))


def main():
    have = os.path.getsize(DEST) if os.path.exists(DEST) else 0
    resume_from = (have // CHUNK) * CHUNK
    if have and have != resume_from:
        with open(DEST, "r+b") as f:
            f.truncate(resume_from)
    todo = []
    pos = resume_from
    while pos < SIZE:
        todo.append((pos, min(pos + CHUNK, SIZE) - 1))
        pos += CHUNK
    print(f"resuming from {resume_from} ({resume_from/1e9:.2f} GB); {len(todo)} chunks", flush=True)
    for rnd in range(1, 4):
        if not todo:
            break
        failed = []
        with ThreadPoolExecutor(max_workers=WORKERS) as ex:
            for start, ok, err in ex.map(job, todo):
                if not ok:
                    failed.append((start, min(start + CHUNK, SIZE) - 1))
                    print(f"  round {rnd} chunk {start} failed: {err}", flush=True)
        print(f"round {rnd}: {len(todo) - len(failed)}/{len(todo)} chunks ok", flush=True)
        todo = failed
    final = os.path.getsize(DEST)
    print(f"DONE size={final} expected={SIZE} match={final == SIZE} failed_chunks={len(todo)}", flush=True)
    if todo:
        print("FAILED offsets:", [t[0] for t in todo], flush=True)


if __name__ == "__main__":
    main()
