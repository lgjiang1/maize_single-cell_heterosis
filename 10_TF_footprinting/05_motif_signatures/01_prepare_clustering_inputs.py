#!/usr/bin/env python
"""
Prepare inputs for motif clustering from JASPAR2026 metadata.

Reads:  5_motif_scanning/plant_motif_family_assignments.tsv  (all 1526 versions)
Writes:
  5_motif_scanning/signatures/motif_names.txt       — one ID per line (latest version per base_id)
  5_motif_scanning/signatures/family_assignments.tsv — family-resolved, one row per selected motif
  5_motif_scanning/signatures/clustering_config.tsv  — per-family motif count (for diagnostics)
"""

import csv
from collections import defaultdict
from pathlib import Path

PROJECT  = Path(__file__).resolve().parents[2]
SCAN_DIR = PROJECT / "5_motif_scanning"
SIG_DIR  = SCAN_DIR / "signatures"
SIG_DIR.mkdir(exist_ok=True)

INPUT  = SCAN_DIR / "plant_motif_family_assignments.tsv"
NAMES  = SIG_DIR / "motif_names.txt"
ASSIGN = SIG_DIR / "family_assignments.tsv"
CONFIG = SIG_DIR / "clustering_config.tsv"

# Fallback: when tf_family is empty, use tf_class (cleaned)
CLASS_TO_FAMILY = {
    "TCP": "TCP",
    "CPP": "CPP",
    "EIL": "EIL",
    "LEAFY": "LEAFY",
    "BBR/BPC": "BBR_BPC",
    "RWP-RK": "RWP-RK",
    "Basic leucine zipper factors (bZIP)": "bZIP",
    "Basic helix-loop-helix factors (bHLH)": "bHLH",
    "C2H2 zinc finger factors": "C2H2",
    "HC3 zinc ribbon factors": "HC3",
    "AP2/EREBP": "AP2",
}

# Normalize inconsistent JASPAR family names
FAMILY_NORMALIZE = {
    "group A": "bZIP_Group_A",
    "D": "bZIP_Group_D",
    "S": "bZIP_Group_S",
    "Group A": "bZIP_Group_A",
    "Group B": "bZIP_Group_B",
    "Group C": "bZIP_Group_C",
    "Group D": "bZIP_Group_D",
    "Group G": "bZIP_Group_G",
    "Group H": "bZIP_Group_H",
    "Group I": "bZIP_Group_I",
    "Group K": "bZIP_Group_K",
    "Group S": "bZIP_Group_S",
    # MADS-box disambiguations
    "Type II": "MADS_Type_II",
    "MIKC": "MADS_MIKC",
}


def main():
    # Read all rows
    rows = []
    with open(INPUT) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rows.append(row)

    # Select latest version per base_id
    latest = {}
    for row in rows:
        bid = row["base_id"]
        ver = int(row["motif_version"].split(".")[-1])
        if bid not in latest or ver > latest[bid][0]:
            latest[bid] = (ver, row)

    # Resolve families and filter
    selected = []
    for bid, (ver, row) in sorted(latest.items()):
        fam = row["tf_family"]
        if not fam:
            cls = row["tf_class"]
            fam = CLASS_TO_FAMILY.get(cls, cls if cls else "")
        if not fam:
            print(f"  DROPPING (no family): {row['motif_version']} {row['tf_name']}")
            continue
        # Normalize and clean family name
        fam = FAMILY_NORMALIZE.get(fam, fam)
        fam_clean = fam.replace("/", "_")
        row["resolved_family"] = fam_clean
        selected.append(row)

    print(f"Selected {len(selected)} motifs (latest version per base_id, family-resolved)")

    # Write motif names list
    with open(NAMES, "w") as f:
        for row in selected:
            # Format: MA0001.2_AGL3 (matches At pipeline convention)
            f.write(f"{row['motif_version']}_{row['tf_name']}\n")
    print(f"Wrote {NAMES}")

    # Write family assignments
    with open(ASSIGN, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["motif_name", "base_id", "motif_version", "tf_name",
                     "tf_family", "tf_class", "species", "tax_id"])
        for row in selected:
            w.writerow([
                f"{row['motif_version']}_{row['tf_name']}",
                row["base_id"],
                row["motif_version"],
                row["tf_name"],
                row["resolved_family"],
                row["tf_class"],
                row["species"],
                row["tax_id"],
            ])
    print(f"Wrote {ASSIGN}")

    # Write per-family config (diagnostics)
    fam_counts = defaultdict(list)
    for row in selected:
        fam_counts[row["resolved_family"]].append(row["tf_name"])

    with open(CONFIG, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["family", "n_motifs", "members"])
        for fam, members in sorted(fam_counts.items(), key=lambda x: -len(x[1])):
            w.writerow([fam, len(members), ";".join(members)])
    print(f"Wrote {CONFIG}")

    # Summary
    print(f"\n{'Family':<30s} {'N':>4s}")
    print("-" * 36)
    for fam, members in sorted(fam_counts.items(), key=lambda x: -len(x[1])):
        print(f"  {fam:<28s} {len(members):4d}")
    print("-" * 36)
    print(f"  {'TOTAL':<28s} {len(selected):4d}")


if __name__ == "__main__":
    main()
