#!/bin/bash

# Input
PATH_VCF_GZ=$1
CHR=$2
PREFIX=$3
DIR_OUTPUT=$4

# Paths
PATH_REF_FASTA=/lizardfs/guarracino/pggb-paper/vcfs/hsapiens/grch38.fa
PATH_REF_SDF=/lizardfs/guarracino/pggb-paper/vcfs/hsapiens/grch38.fa.sdf
DIR_TRUTH_VCF_AND_BED=/lizardfs/guarracino/pggb-paper/vcfs/hsapiens

# Clean up
rm -rf /scratch/$PREFIX
mkdir -p /scratch/$PREFIX
cd /scratch/$PREFIX || exit

# Compare query/truth
for SAMPLE in `zgrep '^#CHROM' $PATH_VCF_GZ -m 1 | cut -f 10- | grep chm13 -v`; do
  echo $SAMPLE
  PATH_TRUTH_VCF_GZ=$DIR_TRUTH_VCF_AND_BED/$SAMPLE.GRCh38_no_alt.deepvariant.pansn.vcf.gz

  TRUTH_VCF_GZ=$(basename $PATH_TRUTH_VCF_GZ .vcf.gz).norm.max50.vcf.gz
  bash /lizardfs/guarracino/pggb-paper/scripts/vcf_preprocess.sh \
    $PATH_TRUTH_VCF_GZ \
    $SAMPLE \
    $PATH_REF_FASTA \
    $CHR \
    50
  # TRUTH_SNV_VCF_GZ=$(basename $TRUTH_VCF_GZ .vcf.gz).snvs.vcf.gz
  # bcftools view --threads 48 -v snps $TRUTH_VCF_GZ -Oz -o $TRUTH_SNV_VCF_GZ && bcftools index -t $TRUTH_SNV_VCF_GZ
  # TRUTH_INDEL_VCF_GZ=$(basename $TRUTH_VCF_GZ .vcf.gz).indels.vcf.gz
  # bcftools view --threads 48 -v indels $TRUTH_VCF_GZ -Oz -o $TRUTH_INDEL_VCF_GZ && bcftools index -t $TRUTH_INDEL_VCF_GZ
  
  QUERY_VCF_GZ=$(basename $PATH_VCF_GZ .vcf.gz).norm.max50.vcf.gz
  bash /lizardfs/guarracino/pggb-paper/scripts/vcf_preprocess.sh \
    $PATH_VCF_GZ \
    $SAMPLE \
    $PATH_REF_FASTA \
    $CHR \
    50
  # QUERY_SNV_VCF_GZ=$(basename $QUERY_VCF_GZ .vcf.gz).snvs.vcf.gz
  # bcftools view --threads 48 -v snps $QUERY_VCF_GZ -Oz -o $QUERY_SNV_VCF_GZ && bcftools index -t $QUERY_SNV_VCF_GZ
  # QUERY_INDEL_VCF_GZ=$(basename $QUERY_VCF_GZ .vcf.gz).indels.vcf.gz
  # bcftools view --threads 48 -v indels $QUERY_VCF_GZ -Oz -o $QUERY_INDEL_VCF_GZ && bcftools index -t $QUERY_INDEL_VCF_GZ
  
  for XXX in "" ".confident-yes" ".easy" ".hard" ".lowmap-segdup"; do
      if [[ "$XXX" == ".confident-yes" ]]; then
          YYY="-e $DIR_TRUTH_VCF_AND_BED/$SAMPLE.f1_assembly_v2.dip$XXX.pansn.bed.gz "
      elif [[ "$XXX" == ".easy" ]]; then
          YYY="-e $DIR_TRUTH_VCF_AND_BED/GRCh38_notinalldifficultregions.bed.gz "
      elif [[ "$XXX" == ".hard" ]]; then
          YYY="-e $DIR_TRUTH_VCF_AND_BED/GRCh38_alldifficultregions.bed.gz "
      elif [[ "$XXX" == ".lowmap-segdup" ]]; then
          YYY="-e $DIR_TRUTH_VCF_AND_BED/GRCh38_alllowmapandsegdupregions.bed.gz "
      else
          YYY=" "
      fi
      SUFFIX="$XXX"
      
      # SNVs+INDELs
      rtg vcfeval \
        -t $PATH_REF_SDF \
        -b $TRUTH_VCF_GZ \
        -c $QUERY_VCF_GZ \
        $YYY-T 48 \
        -o vcfeval/$SAMPLE$SUFFIX

      # # SNVs
      # rtg vcfeval \
      #   -t $PATH_REF_SDF \
      #   -b $TRUTH_SNV_VCF_GZ \
      #   -c $QUERY_SNV_VCF_GZ \
      #   $YYY-T 48 \
      #   -o vcfeval/$SAMPLE.snvs$SUFFIX

      # # INDELs
      # rtg vcfeval \
      #   -t $PATH_REF_SDF \
      #   -b $TRUTH_INDEL_VCF_GZ \
      #   -c $QUERY_INDEL_VCF_GZ \
      #   $YYY-T 48 \
      #   -o vcfeval/$SAMPLE.indels$SUFFIX
  done
done

cd /scratch/
mkdir -p $DIR_OUTPUT
mv /scratch/$PREFIX $DIR_OUTPUT
