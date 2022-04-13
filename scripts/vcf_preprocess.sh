#!/bin/bash

# From https://github.com/wwliao/pangenome-utils/blob/main/preprocess_vcf.sh
# Usage: preprocess_vcf.sh <VCF file> <sample name> <max variant size>

VCF=$1
FNAME=$(basename $VCF)
PREFIX=$(dirname $VCF)/"${FNAME%.vcf.gz}"
SAMPLE=$2
CHROMS=`seq 1 12`
MAXSIZE=$3
REF="SL5.fasta.gz"
MEM="10G"

bcftools=bcftools


#filter
$bcftools view -a -s ${SAMPLE} -Ou ${VCF} \
    | $bcftools norm -f ${REF} -c s -m - -Ou \
    | $bcftools view -e 'GT="ref" | GT~"\."' -f 'PASS,.' -Ou \
    | $bcftools sort -m ${MEM} -T bcftools-sort.XXXXXX -Ou \
    | $bcftools norm -d exact -Oz -o ${PREFIX}.${SAMPLE}.norm.vcf.gz \
    && $bcftools index -t ${PREFIX}.${SAMPLE}.norm.vcf.gz \
    && $bcftools view -e "STRLEN(REF)>${MAXSIZE} | STRLEN(ALT)>${MAXSIZE}" \
                 -r ${CHROMS} -Oz -o ${PREFIX}.max${MAXSIZE}.${SAMPLE}.chr1-12.vcf.gz \
                 ${PREFIX}.${SAMPLE}.norm.vcf.gz \
    #$bcftools index -t ${PREFIX}.max${MAXSIZE}.${SAMPLE}.chr1-12.vcf.gz \
    rm ${PREFIX}.${SAMPLE}.norm.vcf.gz*

## haploid to diploid
    zcat ${PREFIX}.max${MAXSIZE}.${SAMPLE}.chr1-12.vcf.gz|sed "s#GT\t1#GT\t1/1#g" > ${PREFIX}.max${MAXSIZE}.${SAMPLE}.chr1-12.dip.vcf
    bgzip -@ 12 ${PREFIX}.max${MAXSIZE}.${SAMPLE}.chr1-12.dip.vcf
    tabix ${PREFIX}.max${MAXSIZE}.${SAMPLE}.chr1-12.dip.vcf.gz
    rm ${PREFIX}.max${MAXSIZE}.${SAMPLE}.chr1-12.vcf.gz

## snps / indels
$bcftools view -v snps ${PREFIX}.max${MAXSIZE}.${SAMPLE}.chr1-12.dip.vcf.gz -O z -o ${SAMPLE}.pggb.snps.vcf.gz
tabix ${SAMPLE}.pggb.snps.vcf.gz
$bcftools view -v indels ${PREFIX}.max${MAXSIZE}.${SAMPLE}.chr1-12.dip.vcf.gz -O z -o ${SAMPLE}.pggb.indels.vcf.gz
tabix ${SAMPLE}.pggb.indels.vcf.gz
