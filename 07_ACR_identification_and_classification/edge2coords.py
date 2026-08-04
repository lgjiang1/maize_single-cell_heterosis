# This script reads a SAM file containing alignments of ACR edge sequences aligned to the target genome. Each read corresponds to either the left or right edge of an ACR.
# For each alignment that passes the mapping quality (MAPQ) and coverage (QCOV):
# thresholds:
#   - It computes the reference position corresponding to the inner and outer ends of the query sequence (EDGE bp long).
#   - "inner_pos" is the reference coordinate closest to the ACR center.
#   - "outer_pos" is the coordinate farthest away from the ACR center.
#   - It outputs a table with peak ID, side (S or E), target chromosome,inner and outer positions, MAPQ, and query coverage.

# This information will be used in downstream steps to classify whether an unmapped summit corresponds to 
# an insertion, deletion, or more complex structural variation.

# USAGE: python edge2coords.py EDGE MAPQ QCOV < input.sam > output.coords

#!/usr/bin/env python
import sys, re

def die(msg, code=1):
    sys.stderr.write(f"[edge2coords] {msg}\n")
    sys.exit(code)

# ---------------- args ----------------
if len(sys.argv) != 4:
    die("Usage: python edge2coords.py EDGE MAPQ QCOV < in.sam > out.tsv")

try:
    EDGE = int(sys.argv[1])
    MAPQ = int(sys.argv[2])
    QCOV = float(sys.argv[3])
except Exception:
    die("EDGE must be int, MAPQ int, QCOV float")

# -------------- helpers --------------
def qcov_len(cigar: str) -> int:
    cov = 0
    for L, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        L = int(L)
        if op in 'M=XI':
            cov += L
    return cov

def map_qpos_to_ref_fwd(rstart: int, cigar: str, want_qpos: int):
    q = 1
    r = rstart
    for L, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        L = int(L)
        if op in 'M=X': 
            if q + L - 1 >= want_qpos:
                return r + (want_qpos - q)
            q += L; r += L
        elif op == 'I':  
            if q + L - 1 >= want_qpos:
                return r
            q += L
        elif op in 'DN': 
            r += L
        elif op == 'S':  
            if q + L - 1 >= want_qpos:
                return r
            q += L
    return None

def ref_end_from_cigar(rstart: int, cigar: str) -> int:
    """Rightmost 1-based reference coordinate (inclusive)."""
    r = rstart
    for L, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        L = int(L)
        if op in 'M=XDN':
            r += L
    return r - 1

def map_qpos_to_ref_rev(rstart: int, cigar: str, want_qpos: int):
    """
    Map 1-based query position to reference coord accounting for reverse strand.
    We compute total query length (M/= /X, I, and S) then convert want_qpos
    from right-end counting to left-end counting.
    """
    tokens = [(int(L), op) for L, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]
    qlen = sum(L for L, op in tokens if (op in 'M=XI' or op == 'S'))
    target = qlen - want_qpos + 1  

    q = 1
    r = rstart
    for L, op in tokens:
        if op in 'M=X':
            if q + L - 1 >= target:
                return r + (target - q)
            q += L; r += L
        elif op == 'I':
            if q + L - 1 >= target:
                return r
            q += L
        elif op in 'DN':
            r += L
        elif op == 'S':
            if q + L - 1 >= target:
                return r
            q += L
    return None

# =================================================================
# --------------- main -----------------
# =================================================================
print("\t".join(["peakID","side","tchr","pos1","pos2","mapq","qcov"]))

for line in sys.stdin:
    if not line or line[0] == '@':
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 6:
        continue

    qname = f[0]
    try:
        flag  = int(f[1])
        rname = f[2]
        pos   = int(f[3])
        mapq  = int(f[4])
        cigar = f[5]
    except Exception:
        continue

    if rname == '*' or (flag & 0x4):
        continue

    try:
        pid, side_tag = qname.split('|', 1)
    except ValueError:
        continue

    side = 'L' if side_tag.startswith('L') else ('R' if side_tag.startswith('R') else None)
    if side is None:
        continue

    # filters
    cov = qcov_len(cigar) / float(EDGE)
    if mapq < MAPQ or cov < QCOV:
        continue

    is_rev = bool(flag & 16)

    if side == 'L':
        qpos1, qpos2 = 1, EDGE
    else:  # 'R'
        qpos1, qpos2 = 1, EDGE

    if not is_rev:
        pA = map_qpos_to_ref_fwd(pos, cigar, qpos1)
        pB = map_qpos_to_ref_fwd(pos, cigar, qpos2)
    else:
        pA = map_qpos_to_ref_rev(pos, cigar, qpos1)
        pB = map_qpos_to_ref_rev(pos, cigar, qpos2)

    if pA is None or pB is None:
        continue

    if side == 'L':
        outer_ref, inner_ref = pA, pB  # q=1 (outer), q=EDGE (inner)
        if outer_ref > inner_ref:
            outer_ref, inner_ref = inner_ref, outer_ref
        pos1, pos2 = outer_ref, inner_ref
    else:  # 'R'
        inner_ref, outer_ref = pA, pB  # q=1 (inner), q=EDGE (outer)
        if inner_ref > outer_ref:
            inner_ref, outer_ref = outer_ref, inner_ref
        pos1, pos2 = inner_ref, outer_ref

    print("\t".join([pid, side, rname, str(pos1), str(pos2), str(mapq), f"{cov:.3f}"]))
