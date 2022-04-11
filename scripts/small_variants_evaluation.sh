#!/bin/bash


SAMPLE=$1
thread=$2
SDF=$3
# easy region
rtg vcfeval \
        -t $sdf \ 
        -b ../hifi/variant/${SAMPLE}.hifi.snps.vcf.gz \
        -c ${SAMPLE}.pggb.snps.vcf.gz \
        --region=easy_region.bed \
        -T ${thread} \
        -o ${SAMPLE}.snps.easy.report

# hard region
rtg vcfeval \
        -t $SDF \ 
        -b ../hifi/variant/${SAMPLE}.hifi.snps.vcf.gz \
        -c ${SAMPLE}.pggb.snps.vcf.gz \
        --region=hard_region.bed \
        -T ${thread} \
        -o ${SAMPLE}.snps.hard.report
