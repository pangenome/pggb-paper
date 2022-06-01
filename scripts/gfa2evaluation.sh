#!/bin/bash

# This assumes names following the PanSN-spec

PATH_GFA=$1
PREFIX_REFERENCE=$2 #chm13
DIR_OUTPUT=$3
PATH_VCF_PREPROCESS=$4
PATH_NUCMER_2_VCF=$5
THREADS=$6

PREFIX=$(basename "$PATH_GFA" .gfa)

echo "Extracting FASTA file"
PATH_SEQUENCES_FA="$PREFIX".fa
odgi paths -i "$PATH_GFA" -t "$THREADS" -f | bgzip -c -@ "$THREADS" > "$PATH_SEQUENCES_FA".gz
samtools faidx "$PATH_SEQUENCES_FA".gz

echo "--- Take reference sequences"
PATH_REF_FA="$PREFIX"."$PREFIX_REFERENCE".fa
samtools faidx "$PREFIX".fa.gz "$(grep "$PREFIX_REFERENCE" "$PREFIX".fa.gz.fai | cut -f 1)" > "$PATH_REF_FA"

echo "Identify variants with vg"

echo "--- Call variants by haplotype"
# TO FIX: vg deconstruct does not support that mode yet

PATH_SED_GFA="$PREFIX".sed.gfa

# sed replaces only the first instance on a line by default (without the /g modifier)
# To have names like NA21309-1#1#JAHEPC010000450.1 and call haploid genotypes with -H
sed 's/#/-/' "$PATH_GFA" | sed 's/#/#1#/' > "$PATH_SED_GFA"

PATH_VCF="$PREFIX"."$PREFIX_REFERENCE".haplo.vcf
\time -v vg deconstruct -P "$PREFIX_REFERENCE" -H '#' -e -a -t "$THREADS" "$PATH_SED_GFA" > "$PATH_VCF"

# Revert names in the VCF files
grep '^##' "$PATH_VCF" | sed "s/$PREFIX_REFERENCE-1/$PREFIX_REFERENCE/g" > x.vcf
grep '^#CHROM' "$PATH_VCF" | sed 's/-/#/g' >> x.vcf
grep '^#' "$PATH_VCF" -v | sed "s/$PREFIX_REFERENCE-1/$PREFIX_REFERENCE/g" >> x.vcf
mv x.vcf "$PATH_VCF"

echo "--- Take SNVs for each haplotype"
grep '#CHROM' "$PATH_VCF" -m 1 | cut -f 10- | tr '\t' '\n' | while read HAPLO; do
  echo "$HAPLO"

  bash "$PATH_VCF_PREPROCESS" \
      "$PATH_VCF" \
      "$HAPLO" \
      1 \
      "$PATH_REF_FA"
done

echo "Identify variants with nucmer"
NUCMER_VERSION="xxx"
mkdir -p nucmer

echo "--- Align each contig against the reference"
cut -f 1 "$PATH_SEQUENCES_FA".gz.fai | grep "$PREFIX_REFERENCE" -v | while read CONTIG; do
  echo "$CONTIG"

  PREFIX=nucmer/"$CONTIG"
  samtools faidx "$PATH_SEQUENCES_FA".gz "$CONTIG" > "$PREFIX".fa
  echo "$PREFIX" >> tmp
done

cat tmp | parallel -j "$THREADS" "nucmer $PATH_REF_FA {}.fa --prefix {}"
rm tmp

echo "--- Generate VCF files"
cut -f 1 "$PATH_SEQUENCES_FA".gz.fai | grep "$PREFIX_REFERENCE" -v | while read CONTIG; do
  echo "$CONTIG"

  PREFIX=nucmer/"$CONTIG"
  show-snps -THC "$PREFIX".delta > "$PREFIX".var.txt
  #show-snps -TH $PREFIX.delta | cut -f 1-6,9-12 > $PREFIX.var.txt # For taking also variants from ambiguous alignments

  # Check if there are variants
  if [[ $(wc -l "$PREFIX".var.txt | cut -f 1 -d\ ) != 0 ]]; then
    Rscript "$PATH_NUCMER_2_VCF" "$PREFIX".var.txt "$CONTIG" "$PATH_REF_FA" "$NUCMER_VERSION" "$PREFIX".vcf
    bgzip -@ "$THREADS" "$PREFIX".vcf
    tabix "$PREFIX".vcf.gz
  fi
done

echo "--- Merge variants by haplotype"
grep '#CHROM' "$PATH_VCF" -m 1 | cut -f 10- | tr '\t' '\n' | while read HAPLO; do
  echo "$HAPLO"

  grep "^$HAPLO" "$PATH_SEQUENCES_FA".gz.fai | cut -f 1 | while read CONTIG; do
    echo "$CONTIG"

    PREFIX=nucmer/"$CONTIG"
    PATH_VCF=$PREFIX.vcf.gz
    if [[ -f "$PATH_VCF" ]]; then
      # Check if there are valid variants
      if [[ $(zgrep '^#' "$PATH_VCF" -v | grep "^$REF_PREFIX" | wc -l) != 0 ]]; then
        echo "$PATH_VCF" >> "$HAPLO".x.txt
      fi
    fi
  done

  bcftools concat \
    --file-list "$HAPLO".x.txt \
    --naive-force | \
    bcftools reheader --samples <(echo "$HAPLO") > nucmer/"$HAPLO".unsorted.vcf.gz
  rm "$HAPLO".x.txt

  #bcftools sort $HAPLO.unsorted.vcf.gz > $HAPLO.vcf.gz # It doesn't work
  zgrep "^#" nucmer/"$HAPLO".unsorted.vcf.gz > nucmer/"$HAPLO".vcf
  zgrep -v "^#" nucmer/"$HAPLO".unsorted.vcf.gz | sort -k 1,1V -k 2,2n >> nucmer/"$HAPLO".vcf
  bgzip -@ 48 nucmer/"$HAPLO".vcf
  tabix nucmer/"$HAPLO".vcf.gz

  rm nucmer/"$HAPLO".unsorted.vcf.gz
done

echo "Variant evaluation"
echo "--- Prepare the reference in SDF format"
rtg format -o "$PATH_REF_FA".sdf "$PATH_REF_FA"

echo "--- Compare nucmer-based SNVs with PGGB-based SVss"
grep '#CHROM' "$PATH_VCF" -m 1 | cut -f 10- | tr '\t' '\n' | while read HAPLO; do
  echo "$HAPLO"

  PATH_NUCMER_VCF=nucmer/"$HAPLO".vcf.gz
  PATH_PGGB_VCF="$PREFIX"."$PREFIX_REFERENCE".haplo.vcf."$HAPLO".max1.vcf.gz

  # Merge regions closer than 1000 bps to define the callable regions where to evaluate the variants
  dist=1000

  rtg vcfeval \
      -t "$PATH_REF_FA".sdf \
      -b "$PATH_NUCMER_VCF" \
      -c "$PATH_PGGB_VCF" \
      -T "$THREADS" \
      -e <(bedtools intersect -a <(bedtools merge -d $dist -i "$PATH_NUCMER_VCF" ) -b <(bedtools merge -d $dist -i "$PATH_PGGB_VCF")) \
      -o vcfeval/"$HAPLO"
done

cd vcfeval
echo haplotype precision recall f1.score | tr ' ' '\t' > statistics.tsv
grep None */summary.txt | sed 's,/summary.txt:,,' | tr -s ' ' | cut -f 1,7,8,9 -d ' ' | tr ' ' '\t' >> statistics.tsv
cd ..

mkdir -p "$DIR_OUTPUT"
mv "$PREFIX"* nucmer vcfeval "$DIR_OUTPUT"