#!/usr/bin/env python3
"""GSE330794 TRB analysis (pre-registered, replication_design.md).
T1: diversity (Shannon, clonality, top-clone fraction) pre vs post, paired n=7.
T2: clonotype overlap pre/post (Morisita-Horn + shared-clone fractions).
T3: descriptive (P13-post extreme monoclonality, per-sample distributions).
Null: within-patient pre/post label permutation (10k).
"""
import glob, os, re
import numpy as np
import pandas as pd
from scipy import stats

D = "/home/hatch/workspace/bioinfo_downloads/gse330794"
OUT = os.path.join(D, "trb_analysis")
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(7)

def load(f):
    df = pd.read_csv(f, sep="\t")
    # cloneFraction is the authoritative frequency; fall back to counts
    if "cloneFraction" in df.columns:
        p = df["cloneFraction"].values.astype(float)
    else:
        c = df["cloneCount"].values.astype(float)
        p = c / c.sum()
    p = p[p > 0]
    p = p / p.sum()
    return df, p

def shannon(p):
    return float(-(p * np.log(p)).sum())

def clonality(p):
    n = len(p)
    if n <= 1:
        return 1.0
    return float(1.0 - shannon(p) / np.log(n))

def morisita_horn(p, q, keys_p, keys_q):
    # p, q: dicts aaSeqCDR3 -> fraction
    common = set(p) & set(q)
    num = 2.0 * sum(p[k] * q[k] for k in common)
    den = sum(v * v for v in p.values()) + sum(v * v for v in q.values())
    return float(num / den) if den > 0 else 0.0

rows = []
samples = {}
for f in sorted(glob.glob(os.path.join(D, "GSM97329*_TRB.txt.gz"))):
    m = re.search(r"(P\d+)-(pre|post)_TRB", os.path.basename(f))
    pat, tp = m.group(1), m.group(2)
    df, p = load(f)
    n = len(df)
    top = float(p.max())
    top3 = float(np.sort(p)[-3:].sum()) if n >= 3 else float(p.sum())
    rows.append(dict(patient=pat, timepoint=tp, file=os.path.basename(f),
                     n_clones=n, shannon=shannon(p),
                     shannon_norm=shannon(p) / np.log(n) if n > 1 else 0.0,
                     clonality=clonality(p), top_clone=top, top3=top3))
    samples[(pat, tp)] = (df, p)
    print(f"{pat} {tp}: n={n} Shannon={shannon(p):.2f} "
          f"clonality={clonality(p):.3f} top={top:.4f}")

tab = pd.DataFrame(rows)
tab.to_csv(os.path.join(OUT, "trb_diversity.csv"), index=False)

# ---- T1: paired Wilcoxon (n=7) ----
piv = tab.pivot(index="patient", columns="timepoint")
paired = piv.dropna()
print("\npaired patients:", list(paired.index), "n =", len(paired))
for col in ["shannon", "shannon_norm", "clonality", "top_clone", "top3", "n_clones"]:
    post = paired[("shannon", "post")] if False else paired[(col, "post")]
    pre = paired[(col, "pre")]
    d = post - pre
    w = stats.wilcoxon(post, pre)
    print(f"{col}: median post-pre = {d.median():+.4f}, "
          f"up {int((d > 0).sum())}/down {int((d < 0).sum())}, "
          f"Wilcoxon p = {w.pvalue:.4g}")

# permutation null for Shannon (median diff)
post_s = paired[("shannon", "post")].values
pre_s = paired[("shannon", "pre")].values
obs = float(np.median(post_s - pre_s))
nulls = []
for _ in range(10000):
    flip = rng.random(len(post_s)) < 0.5
    d = np.where(flip, pre_s - post_s, post_s - pre_s)
    nulls.append(np.median(d))
nulls = np.array(nulls)
print(f"\nT1 Shannon perm: obs median diff {obs:+.4f}, "
      f"null 95% CI [{np.percentile(nulls, 2.5):+.4f}, {np.percentile(nulls, 97.5):+.4f}], "
      f"empirical p = {(np.abs(nulls) >= abs(obs)).mean():.4f}")

# ---- T2: overlap ----
print("\noverlap (paired patients):")
ov_rows = []
for pat in paired.index:
    df_pre, p_pre = samples[(pat, "pre")]
    df_post, p_post = samples[(pat, "post")]
    d_pre = dict(zip(df_pre["aaSeqCDR3"], p_pre))
    d_post = dict(zip(df_post["aaSeqCDR3"], p_post))
    mh = morisita_horn(d_pre, d_post, None, None)
    # fraction of post repertoire (by abundance) that was present pre
    shared_post = sum(v for k, v in d_post.items() if k in d_pre)
    shared_pre = sum(v for k, v in d_pre.items() if k in d_post)
    n_shared = len(set(d_pre) & set(d_post))
    ov_rows.append(dict(patient=pat, morisita_horn=mh,
                        shared_frac_post=shared_post,
                        shared_frac_pre=shared_pre,
                        n_shared=n_shared,
                        n_pre=len(d_pre), n_post=len(d_post)))
    print(f"{pat}: MH={mh:.4f} shared(post abundance)={shared_post:.3f} "
          f"shared(pre abundance)={shared_pre:.3f} n_shared={n_shared}")
ov = pd.DataFrame(ov_rows)
ov.to_csv(os.path.join(OUT, "trb_overlap.csv"), index=False)
print("\nmedian MH:", float(ov["morisita_horn"].median()))
print("median shared frac of post repertoire:", float(ov["shared_frac_post"].median()))
print("median shared frac of pre repertoire:", float(ov["shared_frac_pre"].median()))

# ---- T3: descriptive ----
print("\nall samples (sorted by clonality):")
print(tab.sort_values("clonality", ascending=False)
      [["patient", "timepoint", "n_clones", "shannon", "clonality", "top_clone"]]
      .to_string(index=False))
print("\nDONE")
