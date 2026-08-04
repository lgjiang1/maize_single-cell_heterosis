"""
07_0a_make_random_regions.py — Per-parent random non-ACR 2kb regions.

For each parent (B73, Ki3, Oh43) sample K=10,000 random 2000bp windows
from chr1-chr10 that do NOT overlap any ACR. These windows serve as the
empirical null for "is this ACR detectably above random?" in step 07_0d.

For B73, the ACR exclusion is the union of B-K and B-O B73-coord ACR sets
(B73 is shared between crosses but each cross has its own peak set).

Reads:
    1_split/genomes/{parent}.fa.fai
    regions/{cross}/{parent}_coords/regions_2000bp.bed   (one per cross)

Writes:
    qc/2_footprinting/random_regions/{parent}_random_2kb.bed   (3-col, no header)

Usage (from project root):
    python 07_0a_make_random_regions.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"

WINDOW_SIZE = 2000
N_WINDOWS = 10000
SEED = 42
KEEP_CHROMS = [f"chr{i}" for i in range(1, 11)]


def load_chrom_sizes(fai_path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with open(fai_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if parts[0] in KEEP_CHROMS:
                sizes[parts[0]] = int(parts[1])
    return sizes


def load_acrs(bed_paths: list[Path]) -> dict[str, np.ndarray]:
    """Return per-chrom sorted (start, end) arrays for the union of all BEDs."""
    frames = []
    for p in bed_paths:
        frames.append(pd.read_csv(
            p, sep="\t", header=None, names=["chrom", "start", "end"]))
    df = pd.concat(frames, ignore_index=True)
    df = df[df["chrom"].isin(KEEP_CHROMS)].copy()

    per_chrom: dict[str, np.ndarray] = {}
    for chrom, grp in df.groupby("chrom"):
        arr = grp[["start", "end"]].to_numpy(dtype=np.int64)
        arr = arr[arr[:, 0].argsort()]
        per_chrom[chrom] = arr
    return per_chrom


def overlaps_any(chrom: str, start: int, end: int,
                 acrs: dict[str, np.ndarray]) -> bool:
    arr = acrs.get(chrom)
    if arr is None or len(arr) == 0:
        return False
    starts, ends = arr[:, 0], arr[:, 1]
    # Candidate ACRs are those with start < end (window's end) and end > start
    # (window's start). Binary search both ends.
    i_lo = np.searchsorted(starts, end, side="left")
    # Among acrs[:i_lo], check if any has end > start
    if i_lo == 0:
        return False
    return bool(np.any(ends[:i_lo] > start))


def sample_for_parent(parent: str, chrom_sizes: dict[str, int],
                      acrs: dict[str, np.ndarray],
                      rng: np.random.Generator) -> pd.DataFrame:
    chroms = list(chrom_sizes.keys())
    lengths = np.array([chrom_sizes[c] for c in chroms], dtype=np.float64)
    # Subtract WINDOW_SIZE so window fits; renormalize.
    avail = np.maximum(lengths - WINDOW_SIZE, 0)
    if avail.sum() == 0:
        raise RuntimeError(f"No chromosome long enough for {WINDOW_SIZE}bp window")
    probs = avail / avail.sum()

    accepted: list[tuple[str, int, int]] = []
    n_tries = 0
    n_rejected = 0
    # Oversample candidates in batches to avoid Python-level loop overhead.
    while len(accepted) < N_WINDOWS:
        batch = max(N_WINDOWS - len(accepted), 1000) * 2
        chosen_chroms = rng.choice(chroms, size=batch, p=probs)
        chosen_starts = rng.integers(
            low=0,
            high=np.array([int(avail[chroms.index(c)]) + 1 for c in chosen_chroms]),
        )
        for c, s in zip(chosen_chroms, chosen_starts):
            n_tries += 1
            s = int(s)
            e = s + WINDOW_SIZE
            if overlaps_any(c, s, e, acrs):
                n_rejected += 1
                continue
            accepted.append((c, s, e))
            if len(accepted) == N_WINDOWS:
                break

    df = pd.DataFrame(accepted, columns=["chrom", "start", "end"])
    # Sort by chrom then start for tidy output.
    df["_chr_order"] = df["chrom"].map({c: i for i, c in enumerate(KEEP_CHROMS)})
    df = (df.sort_values(["_chr_order", "start"])
            .drop(columns="_chr_order")
            .reset_index(drop=True))
    reject_rate = n_rejected / n_tries
    print(f"[INFO] {parent}: {N_WINDOWS} windows kept "
          f"(rejected {n_rejected}/{n_tries} = {100*reject_rate:.1f}% for ACR overlap)")
    print(f"[INFO] {parent}: per-chrom counts:")
    print(df["chrom"].value_counts().reindex(KEEP_CHROMS, fill_value=0).to_string())
    return df


def main() -> int:
    with open(CONFIG_PATH) as f:
        config = yaml.safe_load(f)

    genomes_dir = PROJECT_ROOT / config["paths"]["split_genomes"]
    regions_dir = PROJECT_ROOT / config["paths"]["regions_dir"]
    out_dir = PROJECT_ROOT / "qc" / "2_footprinting" / "random_regions"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Which ACR BEDs feed which parent
    parent_bed_sources: dict[str, list[Path]] = {
        "B73":  [regions_dir / "B-K" / "B73_coords" / "regions_2000bp.bed",
                 regions_dir / "B-O" / "B73_coords" / "regions_2000bp.bed"],
        "Ki3":  [regions_dir / "B-K" / "Ki3_coords" / "regions_2000bp.bed"],
        "Oh43": [regions_dir / "B-O" / "Oh43_coords" / "regions_2000bp.bed"],
    }

    rng = np.random.default_rng(SEED)
    for parent, beds in parent_bed_sources.items():
        for b in beds:
            if not b.exists():
                print(f"[ERROR] missing ACR bed: {b}", file=sys.stderr)
                return 1
        fai = genomes_dir / f"{parent}.fa.fai"
        if not fai.exists():
            print(f"[ERROR] missing fai: {fai}", file=sys.stderr)
            return 1

        print(f"\n[INFO] === {parent} ===")
        chrom_sizes = load_chrom_sizes(fai)
        print(f"[INFO] {parent}: {len(chrom_sizes)} chroms, "
              f"{sum(chrom_sizes.values())/1e6:.0f} Mb total")
        acrs = load_acrs(beds)
        n_acrs = sum(len(v) for v in acrs.values())
        print(f"[INFO] {parent}: excluding {n_acrs:,} ACRs from "
              f"{len(beds)} source BED(s)")

        df = sample_for_parent(parent, chrom_sizes, acrs, rng)
        out_path = out_dir / f"{parent}_random_2kb.bed"
        df.to_csv(out_path, sep="\t", header=False, index=False)
        print(f"[OK] {out_path}")

    print("\n[DONE] random regions written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
