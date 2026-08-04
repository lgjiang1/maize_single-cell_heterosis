#!/usr/bin/env python

import csv, gzip, os
import argparse

#===================================Load sample-to-wells mapping from file============================================
# Supports formats like:
# Pool1 A1-12,B2,B4
# Pool2 E5-7,F1

def load_sample_to_wells(file_path):
    sample_to_wells = {}
    with open(file_path, "r") as f:
        for line in f:
            sample, ranges_str = line.strip().split("\t")
            wells = []
            for r in ranges_str.replace(' ', '').split(','):
                row = r[0]
                body = r[1:]
                if '-' in body:  # eg. A1-12
                    start, end = body.split('-')
                    for i in range(int(start), int(end) + 1):
                        wells.append(f"{row}{i}")
                else:  # eg. A1,F1
                    wells.append(f"{row}{body}")
            sample_to_wells[sample] = wells
    return sample_to_wells

#========================Load barcodes from layout file=====================================================
def load_barcode_layout(layout_file):
    well_to_barcode = {}
    with open(layout_file, "r") as f:
        reader = csv.reader(f, delimiter='\t')
        col_names = next(reader)[1:]
        for row in reader:
            for i, bc in enumerate(row[1:]): 
                well_to_barcode[f"{row[0]}{col_names[i]}"] = bc
    return well_to_barcode

############ Map wells to samples ##########
def build_well_to_sample(sample_to_wells):
    return {w: s for s, wells in sample_to_wells.items() for w in wells}

############ Map barcodes to samples ##########
def build_barcode_to_sample(well_to_barcode, well_to_sample):
    return {bc: well_to_sample[well] for well, bc in well_to_barcode.items() if well in well_to_sample}

#======================================Compute mismatches (N is treated as mismatch)================================
def mismatch_count_with_N(query, ref):
    return sum(0 if q == r else 1 for q, r in zip(query, ref))

#======================================Find best matching barcode (split halves)====================================
#Matches the left and right 5bp separately; both must have ≤1 mismatch
def best_barcode_match_split(query_bc, all_barcodes, max_mismatch=1):
    q_left, q_right = query_bc.split('_')
    best_bc, min_total_mis = None, 999
    for bc in all_barcodes:
        bc_left, bc_right = bc.split('_')
        mis_left = mismatch_count_with_N(q_left, bc_left)
        mis_right = mismatch_count_with_N(q_right, bc_right)
        if mis_left <= max_mismatch and mis_right <= max_mismatch:
            total_mis = mis_left + mis_right
            if total_mis < min_total_mis:
                best_bc, min_total_mis = bc, total_mis  #Keep best match with lowest total mismatches
    return (best_bc, min_total_mis) if best_bc else (None, None)

#========================================Process FASTQ by sample===================================================
def process_fastq(input_fastq, barcode_to_sample, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    all_barcodes = list(barcode_to_sample.keys())
    sample_names = sorted(set(barcode_to_sample.values()))
    base_name = os.path.splitext(os.path.splitext(os.path.basename(input_fastq))[0])[0]
    sample_to_handle = {s: gzip.open(os.path.join(output_dir, f"{base_name}_{s}.fastq.gz"), "wt") for s in sample_names}
    total, n_drop, mis_drop, sample_counts = 0, 0, 0, {s: 0 for s in sample_names}
    with gzip.open(input_fastq, "rt") as fq:
        while (name := fq.readline().rstrip('\n')):
            seq, plus, qual = fq.readline().rstrip('\n'), fq.readline().rstrip('\n'), fq.readline().rstrip('\n')
            total += 1
            bc = "_".join(name.split(' ', 1)[0].split('_')[-2:])
            if bc.count('N') > 2: n_drop += 1; continue  #skip reads with more than 2 'N'
            best, mis = best_barcode_match_split(bc, all_barcodes, max_mismatch=1) #Each half (5bp) is matched separately, allowing up to 1 mismatch per half
            if not best: mis_drop += 1; continue
            s = barcode_to_sample[best]
            sample_to_handle[s].write(f"{name.replace(bc, best, 1)}\n{seq}\n{plus}\n{qual}\n")  #Replace original barcode with corrected one
            sample_counts[s] += 1
    for h in sample_to_handle.values(): h.close()
    #Print summary
    print(f"[{base_name}] total={total}, N={n_drop}, mismatch>{1}={mis_drop}, kept={total - n_drop - mis_drop}")
    for s in sample_names: print(f"  {s}: {sample_counts[s]}")


#=================================================Main function=====================================================
def main():
    parser = argparse.ArgumentParser(description="Demultiplex FASTQ reads by matching split Tn5 barcodes (5+5bp)")
    
    # Define 4 arguments
    parser.add_argument("layout_file", help="96 wells Tn5 barcode layout file path")
    parser.add_argument("input_fastq", help="input clean_reads file path")
    parser.add_argument("sample_well_map", help="Sample-to-well mapping file (e.g., Pool1   A1-12,B1-12)")
    parser.add_argument("--output", default=None, help="Output directory (default: same as input)")
    args = parser.parse_args()

    # Load mappings
    sample_to_wells = load_sample_to_wells(args.sample_well_map)
    well_to_barcode = load_barcode_layout(args.layout_file)
    well_to_sample = build_well_to_sample(sample_to_wells)
    barcode_to_sample = build_barcode_to_sample(well_to_barcode, well_to_sample)

    #Determine output directory
    output_dir = args.output if args.output else os.path.dirname(args.input_fastq)

    # Process the FASTQ
    process_fastq(args.input_fastq, barcode_to_sample, output_dir)


if __name__ == "__main__": 
    main()

