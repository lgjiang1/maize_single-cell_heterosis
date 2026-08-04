#!/usr/bin/bash
#SBATCH --job-name=ACR_classification
#SBATCH --output=logs/classify_ACR_%j.out
#SBATCH --error=logs/classify_ACR_%j.err
#SBATCH --time=10:10:00
#SBATCH --mem=25G
#SBATCH --cpus-per-task=8
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard

########################################
# Configs (Edit names  per run)
########################################
# 1) ref and target genotype name
REF_NAME="B73"
TGT_NAME="Ki3"

# 2) Reference/Target assets
CHAIN="6_ACR_category/B73-Ki3-BxK-KxB/2_cross_genome_analysis/0_input/${REF_NAME}_to_${TGT_NAME}.chain"                     # REF -> TGT
REF_FA="6_ACR_category/B73-Ki3-BxK-KxB/2_cross_genome_analysis/0_input/Zm-${REF_NAME}-REFERENCE-NAM-5.0-chrs-mt-pt.fa"     # REF fasta
REF_GENOME="6_ACR_category/B73-Ki3-BxK-KxB/2_cross_genome_analysis/0_input/${REF_NAME}.genome"                             # REF chrom.sizes
TGT_FA="6_ACR_category/B73-Ki3-BxK-KxB/2_cross_genome_analysis/0_input/Zm-${TGT_NAME}-REFERENCE-NAM-1.0-chrs-mt-pt.fa"     # TGT fasta
PAV_REF="6_ACR_category/B73-Ki3-BxK-KxB/2_cross_genome_analysis/0_input/${REF_NAME}_specific_compared_to_${TGT_NAME}.PAV.bed"  # PAV in REF coords
INPUT_PEAKS="6_ACR_category/B73-Ki3-BxK-KxB/2_cross_genome_analysis/0_input/master_ACR_for_${REF_NAME}-centric.bed"

# 3) Output directory
OUT="6_ACR_category/B73-Ki3-BxK-KxB/2_cross_genome_analysis/1_output"
mkdir -p $OUT

# 4) Peak ID prefix
ID_PREFIX=${REF_NAME}


########################################
# STEP0) generate 1bp summit file
########################################
echo "[0] Make summits for ${REF_NAME}-centric master ACRs"
awk 'BEGIN{OFS="\t"}{print $1,$5-1,$5,"MP"NR,$4,$5}' "$INPUT_PEAKS" > "$OUT/${ID_PREFIX}.summit.bed"


####################################################################
# STEP1) LIFTOVER (Ref SUMMITS -> Target) ==> "Directly-mapped ACR"
####################################################################
echo "LiftOver ${REF_NAME}-centric summits (${REF_NAME} -> ${TGT_NAME})"
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$1,$2,$3}' "$OUT/${ID_PREFIX}.summit.bed" \
| liftOver -bedPlus=3 stdin "$CHAIN" "$OUT/${ID_PREFIX}.summit.on${TGT_NAME}.bed" "$OUT/${ID_PREFIX}.summit.unmapped.bed"

# Write 10-col table (add ACR_type column)
{
  echo -e "${REF_NAME}_chr\t${REF_NAME}_start\t${REF_NAME}_end\t${REF_NAME}_summit\tq_value\t${TGT_NAME}_chr\t${TGT_NAME}_start\t${TGT_NAME}_end\t${TGT_NAME}_summit\tsequence_level"

  awk -v W=250 'BEGIN{OFS="\t"}{
    rchr=$6; rsum=$8; rs=rsum-W-1; if(rs<0) rs=0; re=rsum+W;
    tchr=$1; tsum=$3; ts=tsum-W-1; if(ts<0) ts=0; te=tsum+W;
    print rchr, rs, re, rsum, $5, tchr, ts, te, tsum, "directly-mapped"
  }' "$OUT/${ID_PREFIX}.summit.on${TGT_NAME}.bed"

} > "$OUT/${ID_PREFIX}_directly-mapped-ACR.tsv"


#####################################################################
# STEP2) Identify PAV-ACRs from liftOver-unmapped summits ==> "PAV ACR"
#####################################################################
# Keep REF coords and q-value from the unmapped list
awk 'BEGIN{OFS="\t"} $1!~/^#/ {print $1,$2,$3,$5}' "$OUT/${ID_PREFIX}.summit.unmapped.bed" > "$OUT/${ID_PREFIX}.unmapped.summit.bed"

# Intersect unmapped summits with REF-specific PAV intervals
bedtools intersect -u -a "$OUT/${ID_PREFIX}.unmapped.summit.bed" -b "$PAV_REF" > "$OUT/${ID_PREFIX}.PAV_ACR.summit.bed"

# Make a 10-column table (add ACR_type = "PAV-ACR")
{
  echo -e "${REF_NAME}_chr\t${REF_NAME}_start\t${REF_NAME}_end\t${REF_NAME}_summit\tq_value\t${TGT_NAME}_chr\t${TGT_NAME}_start\t${TGT_NAME}_end\t${TGT_NAME}_summit\tsequence_level"

  awk -v W=250 'BEGIN{OFS="\t"}{
    rchr=$1; rsum=$3; rs=rsum-W-1; if(rs<0) rs=0; re=rsum+W; q=$4;
    print rchr, rs, re, rsum, q, ".", ".", ".", ".", "PAV"
  }' "$OUT/${ID_PREFIX}.PAV_ACR.summit.bed"

} > "$OUT/2_${ID_PREFIX}_PAV-ACR.tsv"


#################################################################################################
# STEP3) Identify SV-rescued ACRs (500bp) from unmapped-not-PAV summits ==> "Structural-shifted ACR"
#################################################################################################
# ------ Parameters ---------
EDGE=150          # Length of flanking edge windows
TOL=10           # Tolerance for insertion/deletion classification on inner span
MAPQ=30           # Minimum mapping quality
QCOV=0.70         # Minimum alignment coverage for edge windows
MAXF=50000       # Upper bound of inner_span to avoid long-distance off-target hits, we think it might be ambiguous insertion
GAP_SMALL_MAX=2000  # Maximum absolute gap size (in bp) to classify as "gap_small", peaks with |delta| ≤ GAP_SMALL_MAX will be kept for downstream DA analysis
WIN=500         # total expected ACR window size

# 3.1 unmapped but not in PAV
bedtools intersect -v -a "$OUT/${ID_PREFIX}.unmapped.summit.bed" -b "$OUT/${ID_PREFIX}.PAV_ACR.summit.bed" > "$OUT/${ID_PREFIX}.unmapped_not_PAV.summit.bed"

# 3.2 expand to ±250 bp and assign PeakIDs
bedtools slop -i "$OUT/${ID_PREFIX}.unmapped_not_PAV.summit.bed" -g "$REF_GENOME" -b 250 \
| awk -v out="$OUT" -v pre="$ID_PREFIX" 'BEGIN{OFS="\t"}{
  id = sprintf(pre "S%06d", NR);
  print $1,$2,$3,id >> out"/"pre".unmapped_not_PAV.peaks_500bp.bed";
  print id,$4        >> out"/"pre".unmapped_not_PAV.id2qvalue.tsv";
}'

PEAKS500="$OUT/${ID_PREFIX}.unmapped_not_PAV.peaks_500bp.bed"

# 3.3 create left/right edge windows (EDGE bp)
awk -v E=$EDGE 'BEGIN{OFS="\t"}{
  chr=$1; s=$2; e=$3; id=$4;
  print chr, s,   s+E, id"|L";
  print chr, e-E, e,   id"|R";
}' "$PEAKS500" > "$OUT/${ID_PREFIX}.unmapped_not_PAV.edges${EDGE}.bed"


# 3.4 extract REF edge sequences
bedtools getfasta -fi "$REF_FA" -bed "$OUT/${ID_PREFIX}.unmapped_not_PAV.edges${EDGE}.bed" -name > "$OUT/${ID_PREFIX}.unmapped_not_PAV.edges${EDGE}.fa"

# 3.5 map edge windows to TGT
minimap2 -x sr -k 15 -w 5 --secondary=no -a -t 8 "$TGT_FA" "$OUT/${ID_PREFIX}.unmapped_not_PAV.edges${EDGE}.fa" > "$OUT/${ID_PREFIX}.unmapped_not_PAV.edges${EDGE}.sam"

# 3.6 parse SAM -> target coords
python ./edge2coords.py "$EDGE" "$MAPQ" "$QCOV" < "$OUT/${ID_PREFIX}.unmapped_not_PAV.edges${EDGE}.sam" > "$OUT/${ID_PREFIX}.unmapped_not_PAV.edges${EDGE}.coords.tsv"

# 3.7 Merge L/R edges to classify and output new intervals
awk 'BEGIN{OFS="\t"}
NR==1{next}   # skip header
{
  k=$1
  side=$2; sub(/::.*/,"",side)    # Replace L::xxx / R::xxx with L / R
  if(side=="L") L[k]=$3"\t"$4"\t"$5
  else if(side=="R") R[k]=$3"\t"$4"\t"$5
  P[k]=1
}
END{
  print "PeakID","L_chr","L_1","L_2","R_chr","R_1","R_2"
  for(k in P){
    split(L[k],a,"\t"); split(R[k],b,"\t")
    print k,a[1],a[2],a[3],b[1],b[2],b[3]
  }
}' "$OUT/${ID_PREFIX}.unmapped_not_PAV.edges${EDGE}.coords.tsv" \
> "$OUT/${ID_PREFIX}.unmapped_not_PAV.edges${EDGE}.coords.wide.tsv"

# 3.8 Classify each peak by edge mapping pattern, identify SV type, and compute the truncated ACR coordinates.
# Output columns:
#   PeakID
#   Category             SV classification label
#   Delta                deviation from expected inner gap (200 bp)
#   Gobs                 observed gaps (R_1 - L_2)   
#   ACR_final_chr        final ACR chromosome
#   ACR_final_start      final ACR start (L_1)
#   ACR_final_end        final ACR end (L_1 + WIN + Delta)
#   ACR_final_len        final ACR length
#   inner_gap_chr        chromosome of the inner gap
#   inner_gap_start      inner gap start (L_2)
#   inner_gap_end        inner gap end (R_1)
#   inner_gap_len        observed inner gap length
#   note                 flags (e.g., L/R missing, diff_chr)

awk -v WIN=$WIN -v EDGE=$EDGE -v TOL=$TOL -v MAXF=$MAXF -v GAP_SMALL_MAX=$GAP_SMALL_MAX -v tgt="$TGT_NAME" '
BEGIN{
  FS=OFS="\t";
  print "PeakID","label","Gobs","delta",
        tgt"_ACR_chr","ACR_final_start","ACR_final_end","ACR_final_len",
        "L_chr","L_1","L_2","R_chr","R_1","R_2","note";
}
function num(x){return (x==""?"":x+0)}
NR==1{next} # skip header
{
  pid=$1; Lchr=$2; L1=num($3); L2=num($4); Rchr=$5; R1=num($6); R2=num($7);
  hasL=(Lchr!=""); hasR=(Rchr!="");
  label=""; delta="NA"; Gobs="NA"; tchr=""; ts=""; te=""; tlen=""; note="";

  # both flanks absent     ## This shouldn’t really happen here in theory, since peaks with both unmapped were already filtered out in Step 6, but I added it just in case :)
  if(!hasL && !hasR){
    label="flank_none"; print pid,label,Gobs,delta,tchr,ts,te,tlen,Lchr,L1,L2,Rchr,R1,R2,"both_unmapped"; next;
  }

  # only left present -> right side missing
  if(hasL && !hasR){
    label="flank_single"; print pid,label,Gobs,delta,tchr,ts,te,tlen,Lchr,L1,L2,Rchr,R1,R2,"R_unmapped"; next;
  }

  # only right present -> left side missing
  if(!hasL && hasR){
    label="flank_single"; print pid,label,Gobs,delta,tchr,ts,te,tlen,Lchr,L1,L2,Rchr,R1,R2,"L_unmapped"; next;
  }

  # cross-chromosome -> ambiguous_translocation
  if(Lchr!=Rchr){
    label="ambiguous_translocation"; print pid,label,Gobs,delta,tchr,ts,te,tlen,Lchr,L1,L2,Rchr,R1,R2,"inter_chr"; next;
  }

  # reverse / inner-cross (R1 < L2) -> inversion
  if(R1!="" && L2!="" && R1<L2){
    Gobs=R1-L2; label="ambiguous_inversion"; tchr=Lchr; ts=L1; te=R2;
    if(te<ts){ tmp=ts; ts=te; te=tmp; note="reverse_span" }
    tlen=te-ts; print pid,label,Gobs,"NA",tchr,ts,te,tlen,Lchr,L1,L2,Rchr,R1,R2,"inner_cross"; next;
  }

  # expected & observed gap
  Gexp=WIN-2*EDGE; Gobs=R1-L2;

  # very large intra-chr gap -> ambiguous_translocation (intra_chr)
  if(Gobs>MAXF){
    label="ambiguous_translocation"; print pid,label,Gobs,"NA",tchr,ts,te,tlen,Lchr,L1,L2,Rchr,R1,R2,"intra_chr"; next;
  }

  delta=Gobs-Gexp; absd=(delta<0?-delta:delta);

  # gap tiering
  if(absd<=TOL) label="minor-shifted";	  
  else if(absd<=GAP_SMALL_MAX) label="small-insertion";
  else label="large-insertion";
  
  # final ACR region (your original definition: start=L1, end=R2)
  tchr=Lchr; ts=L1; te=R2;
  if(te<ts){ tmp=ts; ts=te; te=tmp; label="ambiguous_inversion"; note=(note?note";":"")"outer_reverse" }
  tlen=te-ts;
  print pid,label,Gobs,delta,tchr,ts,te,tlen,Lchr,L1,L2,Rchr,R1,R2,note;
}' "$OUT/${ID_PREFIX}.unmapped_not_PAV.edges${EDGE}.coords.wide.tsv" > "$OUT/${ID_PREFIX}.unmapped_not_PAV.edges${EDGE}.classified.tsv"


# 3.9 Merge BED+qvalue into classified; add missing PeakIDs as flank_none rows
# - File1 (500bp peak BED file) & File2 (qvalue file) & File3 (calssified tsv file)
# - Insert 4 columns (Ref_chr, Ref_start, Ref_end, q_value) right after PeakID.
# - Append missing peaks (present in File1 and File2 but absent in File3) as label=flank_none, other classified columns blank.
gawk -v OFS="\t" -v ref="$REF_NAME" '
BEGIN{FS="\t"}
# 500bp peak BED (id in col4)
ARGIND==1{ if(NF>=4) bed[$4]=$1"\t"$2"\t"$3; next }
# id -> q-value map
ARGIND==2{ if(NF>=2) q[$1]=$2; next }
# classified table (prepend REF coords & q)
ARGIND==3{
  if(FNR==1){
    rest=$0; sub(/^[^\t]+/,"",rest);
    nrest=split(rest,H,"\t");
    print "PeakID",ref"_chr",ref"_start",ref"_end","q_value",rest;
    next
  }
  id=$1; seen[id]=1; split(bed[id],b,"\t");
  chr=(id in bed?b[1]:""); st=(id in bed?b[2]:""); en=(id in bed?b[3]:""); qv=(id in q?q[id]:"");
  printf "%s\t%s\t%s\t%s\t%s",id,chr,st,en,qv;
  for(i=2;i<=NF;i++) printf "\t%s",$i;
  print "";
  next
}
END{
  for(id in bed) if(!(id in seen)){
    split(bed[id],b,"\t"); qv=(id in q?q[id]:"");
    printf "%s\t%s\t%s\t%s\t%s\tflank_none",id,b[1],b[2],b[3],qv;
    for(i=2;i<=nrest;i++) printf "\t";
    print "";
  }
}
' "$OUT/${ID_PREFIX}.unmapped_not_PAV.peaks_500bp.bed" \
  "$OUT/${ID_PREFIX}.unmapped_not_PAV.id2qvalue.tsv" \
  "$OUT/${ID_PREFIX}.unmapped_not_PAV.edges${EDGE}.classified.tsv" \
> "$OUT/${ID_PREFIX}_remaining-ACR.tsv"


####################################################################
# STEP4) Build conserved-ACR set: directlt-mapped + minor-shifted
######################################################################
MINOR_OUT="$OUT/${ID_PREFIX}_minor-shifted-ACR.tsv"

{
  # Write header matching the directly mapped ACR table
  echo -e "${REF_NAME}_chr\t${REF_NAME}_start\t${REF_NAME}_end\t${REF_NAME}_summit\tq_value\t${TGT_NAME}_chr\t${TGT_NAME}_start\t${TGT_NAME}_end\t${TGT_NAME}_summit\tsequence_level"

  awk 'BEGIN{FS=OFS="\t"}
    NR==1{next}
    $6=="minor-shifted"{
      # REF-side coordinates
      ref_chr   = $2
      ref_start = $3
      ref_end   = $4
      q         = $5
      # ACR type
      type      = $6
      # TGT-side coordinates
      tgt_chr   = $9
      tgt_start = $10
      tgt_end   = $11
      print ref_chr, ref_start, ref_end, ".", q,
            tgt_chr, tgt_start, tgt_end, ".", type;
    }
  ' "$OUT/${ID_PREFIX}_remaining-ACR.tsv"

} > "$MINOR_OUT"


# Merge directly-mapped-ACRs and minor-shifted-ACRs into a unified conserved-ACR set
CONSERVE_OUT="$OUT/1_${ID_PREFIX}_Conserved-ACR.tsv"
{
  head -n 1 "$OUT/${ID_PREFIX}_directly-mapped-ACR.tsv"
  {
    tail -n +2 "$OUT/${ID_PREFIX}_directly-mapped-ACR.tsv"
    tail -n +2 "$MINOR_OUT"
  } | sort -t $'\t' -k1,1V -k2,2n
} > "$CONSERVE_OUT"


#############################################################################
# STEP5) Extract insertion-shift ACRs (small-insertion + large-insertion)
#############################################################################
INS_OUT="$OUT/3_${ID_PREFIX}_Insertion-shift-ACR.tsv"

{
  echo -e "${REF_NAME}_chr\t${REF_NAME}_start\t${REF_NAME}_end\t${REF_NAME}_summit\tq_value\t${TGT_NAME}_chr\t${TGT_NAME}_start\t${TGT_NAME}_end\t${TGT_NAME}_summit\tsequence_level"

  awk 'BEGIN{FS=OFS="\t"}
    NR==1{next}   # skip header in remaining-ACR table

    # Keep ACRs labeled as small-insertion or large-insertion
    ($6=="small-insertion" || $6=="large-insertion"){
      ref_chr   = $2    # REF chromosome
      ref_start = $3    # REF start
      ref_end   = $4    # REF end
      q         = $5    # q-value
      type      = $6    # ACR_type (small-insertion / large-insertion)
      tgt_chr   = $9    # TGT chromosome
      tgt_start = $10   # TGT start
      tgt_end   = $11   # TGT end
      # Summit positions are not explicitly defined for these classes;
      # fill both REF_summit and TGT_summit with "."
      print ref_chr, ref_start, ref_end, ".", q,
            tgt_chr, tgt_start, tgt_end, ".", type;
    }
  ' "$OUT/${ID_PREFIX}_remaining-ACR.tsv" \
  | sort -t $'\t' -k1,1V -k2,2n

} > "$INS_OUT"


###############################################################
# STEP6) Extract Ambiguous ACRs (all remaining categories)
###############################################################
AMB_OUT="$OUT/4_${ID_PREFIX}_Ambiguous-ACR.tsv"

{
  head -n 1 "$OUT/${ID_PREFIX}_remaining-ACR.tsv"

  awk 'BEGIN{FS=OFS="\t"}
    NR==1{next} 

    {
      lab = $6
      sub(/\r$/, "", lab)
      gsub(/^[ \t]+|[ \t]+$/, "", lab)
      
      # Keep all rows EXCEPT: minor-shifted, small-insertion, large-insertion
      if (lab!="minor-shifted" && lab!="small-insertion" && lab!="large-insertion") {
        print
      }
    }' "$OUT/${ID_PREFIX}_remaining-ACR.tsv" \
  | sort -t $'\t' -k2,2V -k3,3n

} > "$AMB_OUT"



