#!/usr/bin/env python
"""Build signature_metadata.tsv from Plant_MotifClusters.txt + family_assignments.tsv.

Produces a metadata table with columns:
  signature_id, signature_label, display_name, representative_name, primary_family,
  primary_class, n_members, member_ids, member_names, species_list,
  all_families, is_multi_family

Run after 02a_motif_clustering.R completes.

Usage:
  python 02b_build_signature_metadata.py
"""

import re
from collections import Counter
from pathlib import Path

import pandas as pd


def main():
    project_root = Path(__file__).resolve().parents[2]
    sig_dir = project_root / "5_motif_scanning" / "signatures"

    # Load family assignments
    fam_df = pd.read_csv(sig_dir / "family_assignments.tsv", sep="\t")
    # Build lookup: motif_key (dots→underscores of motif_name) → row
    fam_df["motif_key"] = fam_df["motif_name"].str.replace(".", "_", regex=False)
    fam_lookup = fam_df.set_index("motif_key")
    # Secondary lookup by motif_version key (e.g. MA1253.1 → MA1253_1)
    fam_df["version_key"] = fam_df["motif_version"].str.replace(".", "_", regex=False)
    version_lookup = fam_df.set_index("version_key")

    # Load cluster table
    clusters = {}
    with open(sig_dir / "Plant_MotifClusters.txt") as f:
        for line in f:
            parts = line.strip().split("\t", 1)
            if len(parts) == 2:
                clusters[parts[0]] = parts[1]

    rows = []
    sig_counter = 0
    for sig_label in sorted(clusters.keys(), key=lambda x: int(x.replace("MOTIF", ""))):
        sig_counter += 1
        sig_id = f"sig_{sig_counter:03d}"
        members_str = clusters[sig_label]
        members = [m.strip() for m in members_str.split(";") if m.strip()]

        member_families = []
        member_names = []
        member_ids = []
        member_species = []
        for m in members:
            if m in fam_lookup.index:
                row = fam_lookup.loc[m]
                member_families.append(row["tf_family"])
                member_names.append(row["tf_name"])
                member_ids.append(row["base_id"])
                member_species.append(row["species"])
            else:
                parts = m.split("_")
                version_key = f"{parts[0]}_{parts[1]}" if len(parts) >= 2 else m
                tf_name = "_".join(parts[2:]) if len(parts) >= 3 else m
                if version_key in version_lookup.index:
                    row = version_lookup.loc[version_key]
                    member_families.append(row["tf_family"])
                    member_names.append(row["tf_name"])
                    member_ids.append(row["base_id"])
                    member_species.append(row["species"])
                else:
                    member_names.append(tf_name or m)
                    member_ids.append(m)
                    member_families.append("Unknown")
                    member_species.append("")

        fam_counts = Counter(member_families)
        primary_family = fam_counts.most_common(1)[0][0]
        all_families_set = sorted(set(member_families))
        is_multi_family = len(all_families_set) > 1

        representative = member_names[0] if member_names else "Unknown"

        first_key = members[0] if members[0] in fam_lookup.index else None
        primary_class = ""
        if first_key and first_key in fam_lookup.index:
            pc = fam_lookup.loc[first_key, "tf_class"]
            if pd.notna(pc):
                primary_class = pc

        all_species = sorted(set(s for s in member_species if s))
        display_name = re.sub(r"[/\\: ]+", "_", f"{primary_family}_{representative}")

        rows.append({
            "signature_id": sig_id,
            "signature_label": sig_label,
            "display_name": display_name,
            "representative_name": representative,
            "primary_family": primary_family,
            "primary_class": primary_class,
            "n_members": len(members),
            "member_ids": ";".join(member_ids),
            "member_names": ";".join(member_names),
            "species_list": ";".join(all_species),
            "all_families": ";".join(all_families_set),
            "is_multi_family": is_multi_family,
        })

    df = pd.DataFrame(rows)
    out_path = sig_dir / "signature_metadata.tsv"
    df.to_csv(out_path, sep="\t", index=False)

    print(f"[INFO] Wrote {len(df)} signatures to {out_path}")
    print(f"[INFO] Families: {df['primary_family'].nunique()}")
    print(f"[INFO] Multi-family: {df['is_multi_family'].sum()}")
    print(f"\nFamily distribution:")
    print(df["primary_family"].value_counts().to_string())


if __name__ == "__main__":
    main()
