#!/usr/bin/bash

# ===================== USER PARAMETERS ===================== 

# set two haplotype genome variables
hap1="B73"
hap2="Ki3"

# input ACR tsv files
non_pav="0_input_data/Final_raw_count_for_all_nonPAV-peaks_and_all_samples.tsv"
pav1="0_input_data/Final_raw_count_for_all_PAV-peaks_${hap1}.tsv"
pav2="0_input_data/Final_raw_count_for_all_PAV-peaks_${hap2}.tsv"

# GFF3
gff1="0_input_data/Zm-${hap1}-REFERENCE-NAM-chrs-mt-pt_longest.gff3"
gff2="0_input_data/Zm-${hap2}-REFERENCE-NAM-chrs-mt-pt_longest.gff3"

# Genome sizes
genome1="0_input_data/${hap1}_genome.size"
genome2="0_input_data/${hap2}_genome.size"

# Outputs
outdir="1_filtered_ACR"
mkdir -p $outdir

non_pav_bed="$outdir/Final_raw_count_for_all_nonPAV-peaks_and_all_samples.bedlike.tsv"
out_non_pav="$outdir/Final_raw_count_for_all_nonPAV-peaks_and_all_samples.annot.tsv"
out_pav1="$outdir/Final_raw_count_for_all_PAV-peaks_${hap1}.annot.tsv"
out_pav2="$outdir/Final_raw_count_for_all_PAV-peaks_${hap2}.annot.tsv"

# Thresholds
DIST_CUTOFF=2000   # >2kb intergenic (when no gene overlap)
GENE_FRAC=0.50     # gene_overlap/ACR_len >= 0.5 => genic; otherwise proximal
PART_FRAC=0.50     # exon/intron must cover >0.5 of ACR_len


# ===================== step1: nonPAV --> bedlike =====================
awk -F'\t' -v OFS='\t' '
{
  printf "%s\t%s\t%s\t%s", $2, $3, $4, $1
  for (i = 5; i <= NF; i++) {
    printf "\t%s", $i
  }
  printf "\n"
}
' "$non_pav" > "$non_pav_bed"


# ===================== step2: Define Functions to annotate each ACR based on genomic contex =====================

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  FUNCTION: build features from longest gff3
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
build_features () {
  local gff="$1"
  local genome="$2"
  local featdir="$3"

  mkdir -p "$featdir"

  gawk -F'\t' -v OFS="\t" '
  function attr_get(s,key,  n,i,a,b){
    n=split(s,a,";");
    for(i=1;i<=n;i++){
      split(a[i],b,"=");
      if(b[1]==key) return b[2]
    }
    return ""
  }
  $0 ~ /^#/ {next}

  ($3=="mRNA" || $3=="transcript"){
    tid=attr_get($9,"ID")
    if(tid=="") next
    print $1,$4-1,$5,tid,0,$7 > "'"$featdir"'/tx_body.bed"
    t=($7=="+")?$4:$5
    print $1,t-1,t,tid,0,$7  > "'"$featdir"'/tss.bed"
    next
  }

  $3=="exon"{
    tid=attr_get($9,"Parent")
    if(tid=="") next
    print $1,$4-1,$5,tid > "'"$featdir"'/exon.bed"
  }
  ' "$gff"

  bedtools sort -g "$genome" -i "$featdir/tx_body.bed" > "$featdir/tx_body.sorted.bed"
  bedtools sort -g "$genome" -i "$featdir/exon.bed"    > "$featdir/exon.sorted.bed"
  bedtools sort -g "$genome" -i "$featdir/tss.bed"     > "$featdir/tss.sorted.bed"

  bedtools subtract \
    -a "$featdir/tx_body.sorted.bed" \
    -b "$featdir/exon.sorted.bed" \
  | bedtools sort -g "$genome" -i - \
  > "$featdir/intron.sorted.bed"
}


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FUNCTION: annotate one ACR file
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++
annotate_acr () {
  local acr="$1"
  local genome="$2"
  local featdir="$3"
  local out="$4"

  local tmp
  tmp=$(mktemp -d)
  trap 'rm -rf "$tmp"' RETURN

  head -n 1 "$acr" > "$tmp/header.tsv"
  tail -n +2 "$acr" > "$tmp/body.tsv"

  awk -F'\t' -v OFS="\t" '{print $1,$2-1,$3,$4}' "$tmp/body.tsv" \
  | bedtools sort -g "$genome" -i - \
  > "$tmp/acr.bed"

  bedtools closest -a "$tmp/acr.bed" -b "$featdir/tx_body.sorted.bed" -d -t first -g "$genome" \
  > "$tmp/closest.tsv"

  bedtools intersect -a "$tmp/acr.bed" -b "$featdir/tx_body.sorted.bed" -wo \
  | awk -F'\t' -v OFS="\t" '{s[$4]+=$NF} END{for(i in s) print i,s[i]}' \
  > "$tmp/gene_ov.tsv"

  bedtools intersect -a "$tmp/acr.bed" -b "$featdir/exon.sorted.bed" -wo \
  | awk -F'\t' -v OFS="\t" '{s[$4]+=$NF} END{for(i in s) print i,s[i]}' \
  > "$tmp/exon_ov.tsv"

  bedtools intersect -a "$tmp/acr.bed" -b "$featdir/tss.sorted.bed" -u \
  | awk -F'\t' -v OFS="\t" '{print $4,1}' \
  > "$tmp/tss_hit.tsv"

  awk -F'\t' -v OFS="\t" \
      -v D="$DIST_CUTOFF" \
      -v TG="$GENE_FRAC" \
      -v TH="$PART_FRAC" '
  BEGIN{
    while((getline<"'"$tmp/gene_ov.tsv"'")>0) go[$1]=$2
    while((getline<"'"$tmp/exon_ov.tsv"'")>0) eo[$1]=$2
    while((getline<"'"$tmp/tss_hit.tsv"'")>0) tss[$1]=1
  }
  {
    peak=$4
    len=$3-$2
    dist=$NF

    g=(peak in go)?go[peak]:0
    e=(peak in eo)?eo[peak]:0
    intr=g-e; if(intr<0) intr=0

    if(g==0){
      ann=(dist>D)?"intergenic":"proximal"
    } else if(g/len < TG){
      ann="proximal"
    } else if(tss[peak]){
      ann="TSS"
    } else if(e/len > TH){
      ann="exonic"
    } else if(intr/len > TH){
      ann="intronic"
    } else {
      ann=(e>=intr)?"exonic":"intronic"
    }
    print peak,ann
  }
  ' "$tmp/closest.tsv" \
  | LC_ALL=C sort -t $'\t' -k1,1 \
  > "$tmp/ann.tsv"

  awk -F'\t' -v OFS="\t" '{print $4,$0}' "$tmp/body.tsv" \
  | LC_ALL=C sort -t $'\t' -k1,1 \
  > "$tmp/body.keyed.tsv"

  join -t $'\t' -1 1 -2 1 "$tmp/body.keyed.tsv" "$tmp/ann.tsv" \
  | cut -f2- \
  > "$tmp/body.plus.tsv"

  awk -F'\t' -v OFS="\t" '
    NR==1{
      # header: col1-3 + genomic_annotation + col4..end
      printf "%s\t%s\t%s\tgenomic_annotation", $1,$2,$3
      for(i=4;i<=NF;i++) printf "\t%s", $i
      printf "\n"
      next
    }
    {
      ann=$NF
      # body: col1-3 + ann + col4..(NF-1)
      printf "%s\t%s\t%s\t%s", $1,$2,$3,ann
      for(i=4;i<=NF-1;i++) printf "\t%s", $i
      printf "\n"
    }
  ' "$tmp/header.tsv" "$tmp/body.plus.tsv" > "$out"
}

# ===================== step3: RUN defined functions to generate filtered ACR lists ==============================
echo "Build features"
build_features "$gff1" "$genome1" "$outdir/features_${hap1}"
build_features "$gff2" "$genome2" "$outdir/features_${hap2}"

echo "Annotate ACRs"
annotate_acr "$non_pav_bed" "$genome1" "$outdir/features_${hap1}" "$out_non_pav"
annotate_acr "$pav1" "$genome1" "$outdir/features_${hap1}" "$out_pav1"
annotate_acr "$pav2" "$genome2" "$outdir/features_${hap2}" "$out_pav2"

# To focus specifically on putative cis-regulatory elements, exonic ACRs were excluded
echo "Remove exonic ACRs"
out_non_pav_noex="$outdir/Final_raw_count_for_all_nonPAV-peaks_and_all_samples.annot.noExonic.tsv"
out_pav1_noex="$outdir/Final_raw_count_for_all_PAV-peaks_${hap1}.annot.noExonic.tsv"
out_pav2_noex="$outdir/Final_raw_count_for_all_PAV-peaks_${hap2}.annot.noExonic.tsv"

filter_no_exonic () {
  local in="$1"
  local out="$2"
  awk -F'\t' -v OFS="\t" '
    NR==1 {print; next}
    $4=="exonic" {next}
    {print}
  ' "$in" > "$out"
}

filter_no_exonic "$out_non_pav" "$out_non_pav_noex"
filter_no_exonic "$out_pav1"    "$out_pav1_noex"
filter_no_exonic "$out_pav2"    "$out_pav2_noex"


echo "Removed feature folders"
rm -rf "$outdir/features_${hap1}" "$outdir/features_${hap2}"
rm "$non_pav_bed" "$out_non_pav" "$out_pav1" "$out_pav2"

