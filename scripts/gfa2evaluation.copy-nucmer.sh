#!/bin/bash

# This assumes names following the PanSN-spec

# Example:
#   bash ~/pggb-paper/scripts/gfa2evaluation.sh \
#     /lizardfs/guarracino/day2_hsapiens/graphs/chr20.pan/chr20.pan.fa.gz.254a7b3.04f1c29.5486fb6.smooth.final.gfa \
#     chm13 \
#     out \
#     48

PATH_VCF_PREPROCESS=/home/guarracino/tools/pggb/scripts/vcf_preprocess.sh
PATH_NUCMER_2_VCF=/home/guarracino/tools/pggb/scripts/nucmer2vcf.R

PATH_GFA=$1
PREFIX_REFERENCE=$2
DIR_OUTPUT=$3
THREADS=$4
DIR_NUCMER_TO_COPY=$5



PREFIX=$(basename "$PATH_GFA" .gfa)

echo "Extracting FASTA file"
PATH_SEQUENCES_FA_GZ="$PREFIX".fa.gz
odgi paths -i "$PATH_GFA" -t "$THREADS" -f | bgzip -@ "$THREADS" -l 9  > "$PATH_SEQUENCES_FA_GZ"
samtools faidx "$PATH_SEQUENCES_FA_GZ"

echo "--- Take reference sequences"
PATH_REF_FA="$PREFIX"."$PREFIX_REFERENCE".fa
samtools faidx "$PATH_SEQUENCES_FA_GZ" $(grep "^$PREFIX_REFERENCE""#" "$PATH_SEQUENCES_FA_GZ".fai | cut -f 1) > "$PATH_REF_FA"

echo "Identify variants with vg"

echo "--- Call variants by haplotype"
# TO FIX: vg deconstruct does not support that mode yet

PATH_SED_GFA="$PREFIX".sed.gfa

# sed replaces only the first instance on a line by default (without the /g modifier)
# To have names like NA21309-1#1#JAHEPC010000450.1 and call haploid genotypes with -H
sed 's/#/-/' "$PATH_GFA" | sed 's/#/#1#/' > "$PATH_SED_GFA"

PATH_VCF="$PREFIX"."$PREFIX_REFERENCE".haplo.vcf
PATH_WAVED_VCF="$PREFIX"."$PREFIX_REFERENCE".haplo.waved.vcf
PATH_WAVED_FIXED_VCF="$PREFIX"."$PREFIX_REFERENCE".haplo.waved.fixed.vcf
\time -v vg deconstruct -P "$PREFIX_REFERENCE" -H '#' -e -a -t "$THREADS" "$PATH_SED_GFA" > "$PATH_VCF"
sed -i 's/#0//g' "$PATH_VCF" # Wrong PanSN-spec management in vg 1.41.0

# Revert names in the VCF files
grep '^##' "$PATH_VCF" | sed "s/$PREFIX_REFERENCE-1/$PREFIX_REFERENCE/g" > x.vcf
grep '^#CHROM' "$PATH_VCF" | sed 's/-/#/g' >> x.vcf
grep '^#' "$PATH_VCF" -v | sed "s/^$PREFIX_REFERENCE-1/$PREFIX_REFERENCE/g" >> x.vcf
mv x.vcf "$PATH_VCF"
bgzip -@ "$THREADS" -l 9 "$PATH_VCF"

vcfbub -l 0 -a 10000 --input "$PATH_VCF".gz | vcfwave -I 1000 -t "$THREADS" | bgzip -@ "$THREADS" -l 9 > "$PATH_WAVED_VCF".gz

# The TYPE info sometimes is wrong/missing
# There are variants without the ALT allele
bcftools annotate -x INFO/TYPE "$PATH_WAVED_VCF".gz  | awk '$5 != "."' | bgzip -@ "$THREADS" -l 9 > "$PATH_WAVED_FIXED_VCF".gz


echo "--- Take SNVs for each haplotype"
zgrep '#CHROM' "$PATH_VCF".gz -m 1 | cut -f 10- | tr '\t' '\n' | while read HAPLO; do
  echo "$HAPLO"

  bash "$PATH_VCF_PREPROCESS" \
      "$PATH_VCF".gz \
      "$HAPLO" \
      1 \
      "$PATH_REF_FA"

  bash "$PATH_VCF_PREPROCESS" \
      "$PATH_WAVED_FIXED_VCF".gz \
      "$HAPLO" \
      1 \
      "$PATH_REF_FA"
done

echo "Identify variants with nucmer"
cp -r "$DIR_NUCMER_TO_COPY" .

echo "--- Align each contig against the reference"
echo "--- Generate VCF files"
echo "--- Merge variants by haplotype"

echo "Variant evaluation"
echo "--- Prepare the reference in SDF format"
PATH_REF_SDF="$PREFIX"."$PREFIX_REFERENCE".sdf
rtg format -o "$PATH_REF_SDF" "$PATH_REF_FA"

echo "--- Compare nucmer-based SNVs with PGGB-based SNVs"
zgrep '#CHROM' "$PATH_VCF".gz -m 1 | cut -f 10- | tr '\t' '\n' | while read HAPLO; do
  echo "$HAPLO"

  PATH_NUCMER_VCF=nucmer/"$HAPLO".vcf.gz
  PATH_PGGB_VCF="$PREFIX"."$PREFIX_REFERENCE".haplo."$HAPLO".max1.vcf.gz
  PATH_PGGB_WAVED_VCF="$PREFIX"."$PREFIX_REFERENCE".haplo.waved.fixed."$HAPLO".max1.vcf.gz

  # Merge regions closer than 1000 bps to define the callable regions where to evaluate the variants
  dist=1000

  rtg vcfeval \
      -t "$PATH_REF_SDF" \
      -b "$PATH_NUCMER_VCF" \
      -c "$PATH_PGGB_VCF" \
      -T "$THREADS" \
      -e <(bedtools intersect -a <(bedtools merge -d $dist -i "$PATH_NUCMER_VCF" ) -b <(bedtools merge -d $dist -i "$PATH_PGGB_VCF")) \
      -o vcfeval/haplo/"$HAPLO"

  rtg vcfeval \
      -t "$PATH_REF_SDF" \
      -b "$PATH_NUCMER_VCF" \
      -c "$PATH_PGGB_WAVED_VCF" \
      -T "$THREADS" \
      -e <(bedtools intersect -a <(bedtools merge -d $dist -i "$PATH_NUCMER_VCF" ) -b <(bedtools merge -d $dist -i "$PATH_PGGB_WAVED_VCF")) \
      -o vcfeval/haplo.waved/"$HAPLO"
done

echo "--- Collect statistics"

echo haplotype tp.baseline tp.call fp fn precision recall f1.score nucmer.tot pggb.tot nucmer.ratio pggb.ratio | tr ' ' '\t' > statistics.haplo.tsv
echo haplotype tp.baseline tp.call fp fn precision recall f1.score nucmer.tot pggb.tot nucmer.ratio pggb.ratio | tr ' ' '\t' > statistics.haplo.waved.tsv
zgrep '#CHROM' "$PATH_VCF".gz -m 1 | cut -f 10- | tr '\t' '\n' | while read HAPLO; do
  echo "$HAPLO"

  PATH_NUCMER_VCF=nucmer/"$HAPLO".vcf.gz
  PATH_PGGB_VCF="$PREFIX"."$PREFIX_REFERENCE".haplo."$HAPLO".max1.vcf.gz
  PATH_PGGB_WAVED_VCF="$PREFIX"."$PREFIX_REFERENCE".haplo.waved.fixed."$HAPLO".max1.vcf.gz

  NUM_VARIANTS_NUCMER_TOTAL=$(zgrep '^#' -vc "$PATH_NUCMER_VCF")
  NUM_VARIANTS_PGGB_TOTAL=$(zgrep '^#' -vc "$PATH_PGGB_VCF")
  NUM_VARIANTS_PGGB_WAVED_TOTAL=$(zgrep '^#' -vc "$PATH_PGGB_WAVED_VCF")

  NUM_VARIANTS_NUCMER_TP=$(grep None vcfeval/haplo/"$HAPLO"/summary.txt | sed 's,/summary.txt:,,' | tr -s ' ' | cut -f 3 -d ' ')

  NUM_VARIANTS_PGGB_TP=$(grep None vcfeval/haplo/"$HAPLO"/summary.txt | sed 's,/summary.txt:,,' | tr -s ' ' | cut -f 4 -d ' ')
  FP=$(grep None vcfeval/haplo/"$HAPLO"/summary.txt | sed 's,/summary.txt:,,' | tr -s ' ' | cut -f 5 -d ' ')
  FN=$(grep None vcfeval/haplo/"$HAPLO"/summary.txt | sed 's,/summary.txt:,,' | tr -s ' ' | cut -f 6 -d ' ')
  NUM_VARIANTS_NUCMER_EVALUATED=$(echo "$NUM_VARIANTS_NUCMER_TP + $FN" | bc)
  NUM_VARIANTS_PGGB_EVALUATED=$(echo "$NUM_VARIANTS_PGGB_TP + $FP" | bc)

  NUM_VARIANTS_PGGB_WAVED_TP=$(grep None vcfeval/haplo.waved/"$HAPLO"/summary.txt | sed 's,/summary.txt:,,' | tr -s ' ' | cut -f 4 -d ' ')
  FP_WAVED=$(grep None vcfeval/haplo.waved/"$HAPLO"/summary.txt | sed 's,/summary.txt:,,' | tr -s ' ' | cut -f 5 -d ' ')
  FN_WAVED=$(grep None vcfeval/haplo.waved/"$HAPLO"/summary.txt | sed 's,/summary.txt:,,' | tr -s ' ' | cut -f 6 -d ' ')
  NUM_VARIANTS_NUCMER_WAVED_EVALUATED=$(echo "$NUM_VARIANTS_NUCMER_TP + $FN_WAVED" | bc)
  NUM_VARIANTS_PGGB_WAVED_EVALUATED=$(echo "$NUM_VARIANTS_PGGB_WAVED_TP + $FP_WAVED" | bc)


  # xargs trims whitespaces
  grep None vcfeval/haplo/"$HAPLO"/summary.txt | tr -s ' ' | xargs | cut -f 2,3,4,5,6,7,8,9 -d ' ' | tr ' ' '\t' | \
    awk -v haplo=$HAPLO -v nucmer=$NUM_VARIANTS_NUCMER_TOTAL -v pggb=$NUM_VARIANTS_PGGB_TOTAL -v nucmereval=$NUM_VARIANTS_NUCMER_EVALUATED -v pggbeval=$NUM_VARIANTS_PGGB_EVALUATED -v OFS='\t' '{print(haplo, $0, nucmer, pggb, nucmereval/nucmer, pggbeval/pggb)}' >> statistics.haplo.tsv
  grep None vcfeval/haplo.waved/"$HAPLO"/summary.txt |  tr -s ' ' | xargs | cut -f 2,3,4,5,6,7,8,9 -d ' ' | tr ' ' '\t' | \
      awk -v haplo=$HAPLO -v nucmer=$NUM_VARIANTS_NUCMER_TOTAL -v pggb=$NUM_VARIANTS_PGGB_WAVED_TOTAL -v nucmereval=$NUM_VARIANTS_NUCMER_WAVED_EVALUATED -v pggbeval=$NUM_VARIANTS_PGGB_WAVED_EVALUATED -v OFS='\t' '{print(haplo, $0, nucmer, pggb, nucmereval/nucmer, pggbeval/pggb)}' >> statistics.haplo.waved.tsv
done

mkdir -p "$DIR_OUTPUT"
mv "$PREFIX"* statistics.haplo*.tsv nucmer vcfeval "$DIR_OUTPUT"
