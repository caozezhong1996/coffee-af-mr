# -*- coding: utf-8 -*-
"""Parallel range downloader for ukb-b-5237.vcf.gz (301,395,526 bytes).
Resumable: tracks completed 2MB blocks in a sidecar json."""
import json, os, threading, time
import urllib.request

URL = "https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/PJag1s7QVtMd-HMMySmLTB7Ue-5g2jnqIOHTwx7E8knuFwn8i6xm2R0iyYuSGE72/n/ieup4/b/igd/o/ukb-b-5237/ukb-b-5237.vcf.gz"
BASE = r"C:\Users\64915\Documents\kimi\workspace\integrated_v4\coloc_v42"
OUT = os.path.join(BASE, "ukb-b-5237.vcf.gz")
META = os.path.join(BASE, "vcf_blocks.json")
TOTAL = 301_395_526
BLK = 2 * 1024 * 1024
NBLK = (TOTAL + BLK - 1) // BLK
LOCK = threading.Lock()

if os.path.exists(META):
    done = set(json.load(open(META)))
else:
    done = set()
    with open(OUT, "wb") as f:
        f.truncate(TOTAL)

todo = [i for i in range(NBLK) if i not in done]
print(f"blocks {NBLK}, done {len(done)}, todo {len(todo)}", flush=True)

state = {"finished": 0, "bytes": 0}

def worker():
    while True:
        with LOCK:
            if not todo:
                return
            i = todo.pop(0)
        s, e = i * BLK, min((i + 1) * BLK, TOTAL) - 1
        for attempt in range(4):
            try:
                req = urllib.request.Request(URL, headers={"Range": f"bytes={s}-{e}"})
                with urllib.request.urlopen(req, timeout=120) as r:
                    data = r.read()
                if len(data) != e - s + 1:
                    raise IOError(f"short read {len(data)} != {e-s+1}")
                with LOCK:
                    with open(OUT, "r+b") as f:
                        f.seek(s)
                        f.write(data)
                    done.add(i)
                    json.dump(sorted(done), open(META, "w"))
                    state["finished"] += 1
                    state["bytes"] += len(data)
                break
            except Exception:
                time.sleep(2 + 2 * attempt)
        else:
            with LOCK:
                todo.append(i)  # put back

threads = [threading.Thread(target=worker) for _ in range(24)]
t0 = time.time()
for t in threads:
    t.daemon = True
    t.start()
while any(t.is_alive() for t in threads):
    time.sleep(15)
    el = time.time() - t0
    print(f"progress {len(done)}/{NBLK} blocks, +{state['bytes']/1e6:.0f}MB this run, {el:.0f}s", flush=True)
    if el > 260:
        break
print(f"run end: done {len(done)}/{NBLK}", flush=True)
