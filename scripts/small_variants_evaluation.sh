#!/bin/bash

# Usage: small_variants_evaluation.sh <SAMPLE name> <SDF of ref>

SAMPLE=$1
thread=4
SDF=$2

# snp
# easy region
rtg vcfeval \
        -t $sdf \ 
        -b ${SAMPLE}.hifi.snps.vcf.gz \
        -c ${SAMPLE}.pggb.snps.vcf.gz \
        --region=easy_region.bed \
        -T ${thread} \
        -o ${SAMPLE}.snps.easy.report

# hard region
rtg vcfeval \
        -t $SDF \ 
        -b ${SAMPLE}.hifi.snps.vcf.gz \
        -c ${SAMPLE}.pggb.snps.vcf.gz \
        --region=hard_region.bed \
        -T ${thread} \
        -o ${SAMPLE}.snps.hard.report
