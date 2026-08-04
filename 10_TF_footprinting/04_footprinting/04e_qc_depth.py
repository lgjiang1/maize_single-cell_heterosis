"""
04e_qc_depth.py — Post-footprinting QC: verify outputs + flag low-depth conditions.

Reads:
    _logs/04b_jobs.tsv                 — job table (168 rows)
    3_fragments/{...}/*.frags.tsv.gz   — fragment files (for depth counts)
    4_footprinting/{...}/*/            — scPrinter output directories

Writes:
    qc/2_footprinting/depth_flags.tsv  — 168 rows with depth flags + output verification

Usage (from project root):
    python 04d_qc_depth.py
"""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"


def count_fragments(frag_path: str) -> int | None:
    """Count lines in a gzipped fragment file."""
    p = Path(frag_path)
    if not p.exists():
        return None
    n = 0
    with gzip.open(p, "rt") as f:
        for _ in f:
            n += 1
    return n


def main() -> int:
    with open(CONFIG_PATH) as f:
        config = yaml.safe_load(f)

    fp_cfg = config["footprinting"]
    thresh_ok = int(fp_cfg["depth_flag_ok"])
    thresh_low = int(fp_cfg["depth_flag_low"])

    jobs_path = PROJECT_ROOT / "_logs" / "04b_jobs.tsv"
    if not jobs_path.exists():
        print(f"[ERROR] {jobs_path} not found", file=sys.stderr)
        return 1

    jobs = pd.read_csv(jobs_path, sep="\t")
    print(f"[INFO] loaded {len(jobs)} jobs from {jobs_path}")

    rows = []
    n_missing_frag = 0
    n_missing_fp = 0

    for _, job in jobs.iterrows():
        sample_id = job["sample_id"]
        frag_path = job["frag_out"]
        fp_outdir = Path(job["fp_outdir"])

        # Count fragments
        n_frags = count_fragments(frag_path)
        frag_ok = n_frags is not None

        if not frag_ok:
            print(f"[WARN] missing fragments: {frag_path}")
            n_missing_frag += 1

        # Check scPrinter outputs
        printer_ok = (fp_outdir / "printer.h5ad").exists()
        supp_dir = fp_outdir / "printer_supp"
        tfbs_ok = (supp_dir / "TFBS__ALL.h5ad").exists()
        nucbs_ok = (supp_dir / "NucBS__ALL.h5ad").exists()
        fp_ok = (supp_dir / "FP_zscore__ALL.h5ad").exists()

        if not all([printer_ok, tfbs_ok, nucbs_ok, fp_ok]):
            n_missing_fp += 1
            if not printer_ok:
                print(f"[WARN] missing printer.h5ad: {sample_id}")
            if not tfbs_ok:
                print(f"[WARN] missing TFBS: {sample_id}")
            if not nucbs_ok:
                print(f"[WARN] missing NucBS: {sample_id}")
            if not fp_ok:
                print(f"[WARN] missing FP: {sample_id}")

        # Depth flag (based on fragment count * 2 ≈ read count, or use frags directly)
        # We flag on fragment count since that's what scPrinter sees
        if n_frags is not None:
            approx_reads = n_frags * 2  # each fragment = one read pair
            if approx_reads >= thresh_ok:
                depth_flag = "OK"
            elif approx_reads >= thresh_low:
                depth_flag = "LOW_DEPTH"
            else:
                depth_flag = "VERY_LOW_DEPTH"
        else:
            depth_flag = "MISSING"

        rows.append({
            "cross": job["cross"],
            "coord_system": job["coord_system"],
            "sample_id": sample_id,
            "n_fragments": n_frags if n_frags is not None else 0,
            "approx_reads": (n_frags * 2) if n_frags is not None else 0,
            "depth_flag": depth_flag,
            "frag_ok": frag_ok,
            "printer_ok": printer_ok,
            "tfbs_ok": tfbs_ok,
            "nucbs_ok": nucbs_ok,
            "fp_ok": fp_ok,
        })

    df = pd.DataFrame(rows)
    out_dir = PROJECT_ROOT / "qc" / "2_footprinting"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "depth_flags.tsv"
    df.to_csv(out_path, sep="\t", index=False)

    # Summary
    print(f"\n[INFO] aggregated {len(df)} records")
    print(f"[INFO] missing fragments: {n_missing_frag}")
    print(f"[INFO] missing FP outputs: {n_missing_fp}")
    for flag, count in df["depth_flag"].value_counts().items():
        print(f"[INFO] {flag}: {count}")

    flagged = df[df["depth_flag"].isin(["LOW_DEPTH", "VERY_LOW_DEPTH"])]
    if len(flagged) > 0:
        print(f"\n[WARN] {len(flagged)} conditions flagged for low depth:")
        for _, r in flagged.iterrows():
            print(f"  {r['sample_id']}: {r['approx_reads']:,} reads ({r['depth_flag']})")

    print(f"\n[OK] wrote {out_path} ({len(df)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
