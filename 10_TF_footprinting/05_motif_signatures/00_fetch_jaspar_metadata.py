#!/usr/bin/env python
"""
Fetch TF family/class metadata from JASPAR2026 REST API for all plant motifs.

Input:  5_motif_scanning/motif_inventory.tsv  (extracted from MEME headers)
Output: 5_motif_scanning/plant_motif_family_assignments.tsv

Queries https://jaspar.elixir.no/api/v1/matrix/{id}/ for each unique motif version.
Writes one row per motif version with columns matching the Arabidopsis pipeline format.
"""

import csv
import json
import sys
import time
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError

BASE_URL = "https://jaspar.elixir.no/api/v1/matrix"
PROJECT  = Path(__file__).resolve().parents[2]
SCAN_DIR = PROJECT / "5_motif_scanning"

INVENTORY = SCAN_DIR / "motif_inventory.tsv"
OUTPUT    = SCAN_DIR / "plant_motif_family_assignments.tsv"
CACHE     = SCAN_DIR / "jaspar_api_cache.json"

# Rate limiting
DELAY = 0.3  # seconds between requests


def fetch_matrix(matrix_id: str, retries: int = 3) -> dict:
    """Fetch metadata for a single JASPAR matrix."""
    url = f"{BASE_URL}/{matrix_id}/?format=json"
    for attempt in range(retries):
        try:
            req = Request(url, headers={"Accept": "application/json"})
            with urlopen(req, timeout=30) as resp:
                return json.loads(resp.read().decode())
        except (HTTPError, URLError, TimeoutError) as e:
            if attempt < retries - 1:
                wait = 2 ** attempt
                print(f"  retry {attempt+1} for {matrix_id} after {wait}s: {e}")
                time.sleep(wait)
            else:
                print(f"  FAILED {matrix_id}: {e}")
                return None


def main():
    # Read inventory
    motifs = []
    with open(INVENTORY) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                motifs.append({
                    "file_stem": parts[0],      # e.g. MA0001.1
                    "motif_version": parts[1],   # e.g. MA0001.1
                    "tf_name": parts[2],         # e.g. AGL3
                })

    print(f"Loaded {len(motifs)} motif entries from inventory")

    # Load cache if exists
    cache = {}
    if CACHE.exists():
        with open(CACHE) as f:
            cache = json.load(f)
        print(f"Loaded {len(cache)} cached entries")

    # Fetch metadata for each unique motif_version
    unique_versions = sorted(set(m["motif_version"] for m in motifs))
    to_fetch = [v for v in unique_versions if v not in cache]
    print(f"Need to fetch {len(to_fetch)}/{len(unique_versions)} motifs from API")

    for i, mid in enumerate(to_fetch):
        data = fetch_matrix(mid)
        if data:
            cache[mid] = {
                "family": data.get("family", []),
                "class": data.get("class", []),
                "species": [s.get("name", "") for s in data.get("species", [])],
                "tax_id": [s.get("tax_id", "") for s in data.get("species", [])],
            }
        else:
            cache[mid] = {"family": [], "class": [], "species": [], "tax_id": []}

        if (i + 1) % 50 == 0 or i == len(to_fetch) - 1:
            print(f"  fetched {i+1}/{len(to_fetch)}")
            # Save cache periodically
            with open(CACHE, "w") as f:
                json.dump(cache, f, indent=1)

        time.sleep(DELAY)

    # Final cache save
    with open(CACHE, "w") as f:
        json.dump(cache, f, indent=1)
    print(f"Cache saved: {len(cache)} entries")

    # Write output TSV
    with open(OUTPUT, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "motif_name", "base_id", "motif_version", "tf_name",
            "tf_family", "tf_class", "species", "tax_id"
        ])
        for m in motifs:
            ver = m["motif_version"]
            base_id = ver.rsplit(".", 1)[0]
            meta = cache.get(ver, {})
            w.writerow([
                f"{ver}_{m['tf_name']}",        # motif_name (matches At format)
                base_id,                          # base_id
                ver,                              # motif_version
                m["tf_name"],                     # tf_name
                ";".join(meta.get("family", [])), # tf_family
                ";".join(meta.get("class", [])),  # tf_class
                ";".join(meta.get("species", [])),# species
                ";".join(str(x) for x in meta.get("tax_id", [])),  # tax_id
            ])

    print(f"Wrote {len(motifs)} rows to {OUTPUT}")

    # Summary
    families = {}
    for m in motifs:
        ver = m["motif_version"]
        meta = cache.get(ver, {})
        fam = ";".join(meta.get("family", [])) or "UNKNOWN"
        families[fam] = families.get(fam, 0) + 1

    print(f"\nFamily distribution ({len(families)} families):")
    for fam, n in sorted(families.items(), key=lambda x: -x[1])[:25]:
        print(f"  {fam:30s} {n:4d}")


if __name__ == "__main__":
    main()
