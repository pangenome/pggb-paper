#!/bin/bash

# Input
PATH_VCF_GZ=$1
PREFIX=$2
DIR_OUTPUT=$3

# Paths
PATH_REF_FASTA=/lizardfs/guarracino/pggb-paper/vcfs/hsapiens/chm13v2.fa
PATH_REF_SDF=/lizardfs/guarracino/pggb-paper/vcfs/hsapiens/chm13v2.fa.sdf
DIR_TRUTH_VCF=/lizardfs/guarracino/pggb-paper/vcfs/hsapiens

CHR=chm13#1#chr6

# Clean up
rm -rf /scratch/$PREFIX
mkdir -p /scratch/$PREFIX
cd /scratch/$PREFIX || exit

# Compare query/truth
for SAMPLE in `zgrep '^#CHROM' $PATH_VCF_GZ -m 1 | cut -f 10-`; do
  echo $SAMPLE
  PATH_TRUTH_VCF_GZ=$DIR_TRUTH_VCF/$SAMPLE.CHM13Y.deepvariant.vcf.gz

  bash /lizardfs/guarracino/pggb-paper/scripts/vcf_preprocess.sh \
    $PATH_TRUTH_VCF_GZ \
    $SAMPLE \
    $PATH_REF_FASTA \
    $CHR \
    50
  TRUTH_VCF_GZ=$(basename $PATH_TRUTH_VCF_GZ .vcf.gz).norm.max50.vcf.gz

  bash /lizardfs/guarracino/pggb-paper/scripts/vcf_preprocess.sh \
    $PATH_VCF_GZ \
    $SAMPLE \
    $PATH_REF_FASTA \
    $CHR \
    50
  QUERY_VCF_GZ=$(basename $PATH_VCF_GZ .vcf.gz).norm.max50.vcf.gz

  rtg vcfeval \
    -t $PATH_REF_SDF \
    -b $TRUTH_VCF_GZ \
    -c $QUERY_VCF_GZ \
    -T 48 \
    -o vcfeval/$SAMPLE
done

cd /scratch/
mkdir -p $DIR_OUTPUT
mv /scratch/$PREFIX $DIR_OUTPUT
