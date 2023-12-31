#!/bin/bash

# Input
PATH_VCF_GZ=$1            # Path to bgzipped VCF file
CHR=$2                    # Chromosome
PREFIX=$3                 # Prefix for output directory
DIR_OUTPUT=$4             # Directory for output
PATH_REF_FASTA=$5         # Path to reference FASTA
PATH_REF_SDF=$6           # Path to reference SDF
DIR_TRUTH_VCF_AND_BED=$7  # Directory containing truth VCF and BED files
TRUTH_VCF_PATTERN=$8      # Pattern for truth VCF filenames
CONDITIONS=("${@:9}")     # Additional conditions (e.g. ".confident-yes:PATH_TO_BED_GZ")

# Clean up
rm -rf /scratch/$PREFIX
mkdir -p /scratch/$PREFIX
cd /scratch/$PREFIX || exit

# Compare query/truth
for SAMPLE in `zgrep '^#CHROM' $PATH_VCF_GZ -m 1 | cut -f 10- | tr '\t' '\n' | grep chm13 -v | tr '\n' ' '`; do
  echo $SAMPLE
  PATH_TRUTH_VCF_GZ=$DIR_TRUTH_VCF_AND_BED/${SAMPLE}$TRUTH_VCF_PATTERN
  # Check if PATH_TRUTH_VCF_GZ exists
  if [[ ! -f "$PATH_TRUTH_VCF_GZ" ]]; then
      echo "Warning: File $PATH_TRUTH_VCF_GZ not found. Skipping $SAMPLE."
      continue
  fi
  
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
  
  # Loop through each condition
  for CONDITION in "${CONDITIONS[@]}"; do
      # Split CONDITION into SUFFIX and BED_PATH
      SUFFIX="${CONDITION%%:*}"
      BED_PATH="${CONDITION#*:}"
      
      YYY=" "
      if [[ -n "$BED_PATH" ]]; then
          YYY="--bed-regions $BED_PATH "
      fi
      
      echo "$SUFFIX $YYY"
      
      # SNVs+INDELs
      rtg vcfeval \
        -t $PATH_REF_SDF \
        -b $TRUTH_VCF_GZ \
        -c $QUERY_VCF_GZ \
        $YYY-T 48 \
        -o vcfeval/$SAMPLE$SUFFIX \
        -e <(bedtools intersect -a <(bedtools merge -d 1000 -i $TRUTH_VCF_GZ) -b <(bedtools merge -d 1000 -i $QUERY_VCF_GZ))

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
echo "mv /scratch/$PREFIX $DIR_OUTPUT"
mv /scratch/$PREFIX $DIR_OUTPUT
