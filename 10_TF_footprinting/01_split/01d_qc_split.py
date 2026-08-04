"""
01d_qc_split.py — aggregate per-BAM stats from 01b into QC tables.

Walks `_logs/01b_jobs.tsv` and reads each output BAM's `*.stats.json` (written
by 01b_split_and_rename.py), then emits three TSVs under `qc/1_split/`:

  - depth_table.tsv     : per-output-BAM read counts and quality stats
  - crossmap_qc.tsv     : F1 sub-BAM cross-mapping QC (q30 loss fraction)
  - parental_floor.tsv  : per (cross, cell_type, coord_system) parental depth,
                          which Phase 2 will use as the asymmetric subsampling
                          target for the F1 sub-BAMs

Usage (run from project root, after 01b array job completes):
    python 01d_qc_split.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"


def load_config() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def load_jobs(jobs_tsv: Path) -> list[dict]:
    rows: list[dict] = []
    with open(jobs_tsv) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                rows.append({
                    "input_bam": parts[0],
                    "target_parent": parts[1],
                    "output_bam": parts[2],
                })
    return rows


def parse_input_meta(input_bam: str) -> dict:
    """Derive cross / role / name / cell_type from an input BAM path.

    e.g. '0_bams/0_B-K_input/B73-C1.renamed.bam'
        → cross=B-K, role=P1, name=B73, cell_type=C1
    """
    p = Path(input_bam)
    parent_dir = p.parent.name              # '0_B-K_input'
    cross = parent_dir.replace("0_", "").replace("_input", "")  # 'B-K'

    fname = p.name.replace(".renamed.bam", "").replace(".bam", "")
    if "-C" in fname:
        name, cell_num = fname.rsplit("-C", 1)
        cell_type = f"C{cell_num}"
    else:
        name = fname
        cell_type = "?"

    if name == "B73":
        role = "P1"
    elif name in ("Ki3", "Oh43"):
        role = "P2"
    elif name in ("BxK", "BxO"):
        role = "F1_fwd"
    elif name in ("KxB", "OxB"):
        role = "F1_rev"
    else:
        role = "?"

    return {"cross": cross, "role": role, "name": name, "cell_type": cell_type}


def parse_output_meta(output_bam: str) -> dict:
    """Derive coord_system and allele suffix from an output BAM path."""
    p = Path(output_bam)
    coord_system = p.parent.name.replace("_coords", "")  # 'B73_coords' → 'B73'

    fname = p.name.replace(".bam", "")
    allele = ""
    for tag in ("B73allele", "Ki3allele", "Oh43allele"):
        if fname.endswith(f".{tag}"):
            allele = tag
            break
    return {"coord_system": coord_system, "allele": allele}


def main() -> int:
    config = load_config()
    qc_dir = PROJECT_ROOT / config["paths"]["qc_dir"] / "1_split"
    qc_dir.mkdir(parents=True, exist_ok=True)

    jobs_tsv = PROJECT_ROOT / "_logs" / "01b_jobs.tsv"
    if not jobs_tsv.exists():
        print(f"[ERROR] jobs TSV not found: {jobs_tsv}", file=sys.stderr)
        print(f"[ERROR] run 01b_make_jobs.py first", file=sys.stderr)
        return 1

    jobs = load_jobs(jobs_tsv)
    print(f"[INFO] loaded {len(jobs)} jobs from {jobs_tsv.relative_to(PROJECT_ROOT)}")

    # Aggregate per-BAM stats
    records: list[dict] = []
    missing = 0
    for job in jobs:
        out_path = PROJECT_ROOT / job["output_bam"]
        stats_path = out_path.with_suffix(out_path.suffix + ".stats.json")
        if not stats_path.exists():
            missing += 1
            print(f"[WARN] missing stats: {stats_path.relative_to(PROJECT_ROOT)}", file=sys.stderr)
            continue
        with open(stats_path) as f:
            stats = json.load(f)
        in_meta = parse_input_meta(job["input_bam"])
        out_meta = parse_output_meta(job["output_bam"])
        rec = {**in_meta, **out_meta, **stats}
        records.append(rec)

    print(f"[INFO] aggregated {len(records)} records  (missing {missing})")
    if missing > 0:
        print(f"[WARN] {missing}/{len(jobs)} stats files missing — partial output", file=sys.stderr)

    # ----- depth_table.tsv -----
    depth_columns = [
        "cross", "coord_system", "name", "role", "cell_type", "allele",
        "total_reads", "mapped_reads", "primary_reads", "target_parent_reads",
        "mate_ok_reads", "mapq10_pass", "mapq30_pass", "q30_loss_fraction",
        "output_reads", "output_bam",
    ]
    depth_path = qc_dir / "depth_table.tsv"

    def sort_key(r):
        # group by cross, cell type number, coord system, role
        try:
            cell_n = int(r["cell_type"][1:])
        except (ValueError, KeyError):
            cell_n = 999
        return (r["cross"], cell_n, r["coord_system"], r.get("role", ""), r.get("allele", ""))

    sorted_records = sorted(records, key=sort_key)

    with open(depth_path, "w") as f:
        f.write("\t".join(depth_columns) + "\n")
        for rec in sorted_records:
            row = []
            for c in depth_columns:
                v = rec.get(c, "")
                if v is None:
                    v = "NA"
                row.append(str(v))
            f.write("\t".join(row) + "\n")
    print(f"[OK]   wrote {depth_path.relative_to(PROJECT_ROOT)} ({len(records)} rows)")

    # ----- crossmap_qc.tsv (F1 only) -----
    qc_max = float(config["split"]["qc_max_q30_loss"])
    crossmap_columns = [
        "cross", "name", "cell_type", "coord_system", "allele",
        "target_parent_reads", "mate_ok_reads", "mapq10_pass", "mapq30_pass",
        "q30_loss_fraction", "flag",
    ]
    crossmap_path = qc_dir / "crossmap_qc.tsv"
    f1_records = [r for r in sorted_records if r["role"] in ("F1_fwd", "F1_rev")]

    with open(crossmap_path, "w") as f:
        f.write("\t".join(crossmap_columns) + "\n")
        for rec in f1_records:
            loss = rec.get("q30_loss_fraction")
            flag = "OVER_THRESHOLD" if (loss is not None and loss > qc_max) else "ok"
            loss_str = f"{loss:.4f}" if loss is not None else "NA"
            row = [
                rec["cross"], rec["name"], rec["cell_type"], rec["coord_system"], rec.get("allele", ""),
                str(rec["target_parent_reads"]), str(rec["mate_ok_reads"]),
                str(rec["mapq10_pass"]), str(rec["mapq30_pass"]),
                loss_str, flag,
            ]
            f.write("\t".join(row) + "\n")
    print(f"[OK]   wrote {crossmap_path.relative_to(PROJECT_ROOT)} ({len(f1_records)} rows)")

    # ----- parental_floor.tsv -----
    # One row per (cross, cell_type, coord_system) — the parental BAM's output
    # read count is the floor that Phase 2 will subsample F1 sub-BAMs down to.
    floor_columns = ["cross", "cell_type", "coord_system", "parent_name", "parent_reads"]
    floor_path = qc_dir / "parental_floor.tsv"
    parental_records = [r for r in sorted_records if r["role"] in ("P1", "P2")]

    with open(floor_path, "w") as f:
        f.write("\t".join(floor_columns) + "\n")
        for rec in parental_records:
            row = [
                rec["cross"], rec["cell_type"], rec["coord_system"],
                rec["name"], str(rec["mapq10_pass"]),
            ]
            f.write("\t".join(row) + "\n")
    print(f"[OK]   wrote {floor_path.relative_to(PROJECT_ROOT)} ({len(parental_records)} rows)")

    # Surface flagged entries
    flagged = [r for r in f1_records if r.get("q30_loss_fraction") is not None and r["q30_loss_fraction"] > qc_max]
    if flagged:
        print(f"\n[WARN] {len(flagged)} F1 sub-BAMs exceed q30 loss threshold ({qc_max:.0%}):", file=sys.stderr)
        for r in flagged[:20]:
            print(
                f"  {r['cross']:5s} {r['name']:5s} {r['cell_type']:4s} {r['coord_system']:5s}  "
                f"q30_loss={r['q30_loss_fraction']:.4f}",
                file=sys.stderr,
            )
        if len(flagged) > 20:
            print(f"  ... and {len(flagged) - 20} more", file=sys.stderr)
    else:
        print(f"\n[OK]  no F1 sub-BAMs exceed the q30 loss threshold ({qc_max:.0%})")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
