# Soy Pangenome

## Software versions

```
rtg-tools: 3.12.1
mash: 2.2
bcftools: 1.13
samtools: 1.12
```

## Data collection

Download the assemblies from the two studies:
- Liu et al. Pan-Genome of Wild and Cultivated Soybeans. Cell. 2020 Jul 9;182(1):162-176.e13. doi: [10.1016/j.cell.2020.05.023](https://www.sciencedirect.com/science/article/pii/S0092867420306188)
- Chu et al. Eight soybean reference genome resources from varying latitudes and agronomic traits. Sci Data. 2021;8,164. doi: [10.1038/s41597-021-00947-2](https://doi.org/10.1038/s41597-021-00947-2)

```
mkdir assemblies/soy
cd assemblies/soy
cat ../data/soy.urls| parallel -j 4 'wget -q {} && echo got {}'
```

## Pangenome Sequence Naming and Partitioning into Chromosomes

Extract the separate chromosomes from the assemblies, rename them, and discard the unplaced small contigs.

```
for i in WM82 ZH13 C01 C02 C03 C04 C05 C06 C07 C08 C09 C10 C11 C12 C13 C14 L01 L02 L03 L04 L05 L06 L07 L08 L09 W01 W02 W03 W05 GmWF7 GmHF25 GmZH35 GmZH13 GmJY GmHX3 GmW82 GsojaF; do 
    for file in *$i*a.gz; do 
        for c in {01..20}; do
            start=`zgrep -n "^>" $file | head -n$c | tail -n1 | cut -d: -f1`
            end=`zgrep -n "^>" $file | head -n$((c+1)) | tail -n1 | cut -d: -f1`
            zcat $file | \
                head -n $((end-1)) | \
                tail -n $((end-start)) | \
                sed "s/^>.*$/>$i#1#chr$c/" \
                >> SOY37.chr$c.fa
        done
    done
    echo "Processed $i"
done
for c in {01..20}; do
    bgzip -@10 SOY37.chr$c.fa
    samtools faidx SOY37_chr$c.fa.gz
done
```

For efficient testing, we continue only with building a chromosome 18 graph.

## Divergence estimation

Calculate mash distances between the chromosomes to determine a suitable minimum percent sequence identity for the alignments in the graph generation:
```
mash triangle SOY37.chr18.fa | tr " " "\t" > SOY37.chr18.fa.mash_triangle
echo "max distances:"
sed 1,1d SOY37.chr18.fa.mash_triangle | tr '\t' '\n' | grep -v chr | sort -g -k 1nr | head -n 5
```

Sequences are closely related, max mash distance ~ 0.01. Therefore, the pggb parameter 'min. percent identity' will be tested with 0.98 and 0.95. We are a bit careful, since a MDS plot based on euclidean distances (`odgi paths -d`) shows the sample W05 as an extreme outlier, while the mash distances do not.


## Pangenome graph building

```
# used: ghcr.io/pangenome/pggb:latest on May 8th
bgzip SOY37.chr18.fa
samtools faidx SOY37.chr18.fa.gz
# started on a kubernetes cluster via Rancher, but command should be sth like
docker run -it -v /path/to/graph:/data ghcr.io/pangenome/pggb:latest "pggb -i /data/SOY37.chr18.fa.gz -s50000 -p98 -n36 -H37 -k49 -t30 -T25 -G5000,10000 -S -M -o /data/chr18_s50p98k49 -V \"ZH13:#\""
docker run -it -v /path/to/graph:/data ghcr.io/pangenome/pggb:latest "pggb -i /data/SOY37.chr18.fa.gz -s30000 -p95 -n36 -H37 -k47 -t30 -T25 -G5000,10000 -S -M -o /data/chr18_s30p95k47 -V \"ZH13:#\""
```

## Graph statistics
...


## Small variants evaluation

Download variants from the Liu 2020 publication from FigShare (https://figshare.com/s/689ae685ad2c368f2568):
```
wget https://figshare.com/ndownloader/files/22137132\?private_link\=689ae685ad2c368f2568
mv 22137132\?private_link\=689ae685ad2c368f2568 29.pan-genome.SNP.txt

wget https://figshare.com/ndownloader/files/22137192\?private_link\=689ae685ad2c368f2568
mv 22137192\?private_link\=689ae685ad2c368f2568 29.pan-genome.Small-INDEL.txt
```

Convert variants in text files from the Liu 2020 publication into a fake vcf:
```
echo "##fileformat=VCFv4.2" > header
echo "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >> header
echo "##contig=<ID=ZH13#1#chr18,length=60669500>" >> header
cat header > Liu.SNPs.chr18.vcf
awk '$1=="Chr"{for (i=5; i<=NF; ++i) { s=s $i; if (i<NF) s=s "\t" } print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" s} $1=="Chr18"{$1="ZH13#1#chr18"; s=""; for (i=5; i<=NF; ++i) { if (i>5) s=s "\t"; if ($i == ".") s=s "./."; else s=s "1/1"; } print $1 "\t" $2 "\t.\t" $3 "\t" $4 "\t.\tPASS\t.\tGT\t" s}' 29.pan-genome.SNP.txt | sed "s/Wm82/WM82/" | sed "s/Soy//g" >> Liu.SNPs.chr18.vcf

cat header > Liu.INDELs.chr18.vcf
awk '$1=="Chr"{for (i=6; i<=NF; ++i) { s=s $i; if (i<NF) s=s "\t" } print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" s} $1=="Chr18"{$1="ZH13#1#chr18"; s=""; for (i=6; i<=NF; ++i) { if (i>6) s=s "\t"; if ($i == ".") s=s "./."; else s=s "1/1"; } print $1 "\t" $2 "\t.\t" $4 "\t" $5 "\t.\tPASS\t.\tGT\t" s}' 29.pan-genome.Small-INDEL.txt | sed "s/Wm82/WM82/" | sed "s/Soy//g" >> Liu.INDELs.chr18.vcf
```
!!! TODO: multiallelic sites in raw variant files are not considered, typically the major allele is reported as variant allele, the minor alleles discarded. Will have to be fixed for final publication.


Extract variants from the pggb graphs:
```
# for now only for one sample C01
SAMPLE=C01

for PREFIX in chr18_s30p95k47/SOY37.chr18.fa.gz.8a73408.e34d4cd.8bfd93c.smooth.final chr18_s50p98k49/SOY37.chr18.fa.gz.f330e33.f85392c.87069a2.smooth.final; do
    # Andrea style:
    bash pggb/scripts/vcf_preprocess.sh \
        $PREFIX.ZH13.vcf \
        $SAMPLE \       # extract variants from that sample
        1 \             # extract only SNPs
        ZH13.chr18.fa   # reference

    # Joerg style:
    bcftools norm \
        -m -any \
        $PREFIX.ZH13.vcf \
        | bcftools view \
            -s $SAMPLE \
            -v snps \
            --min-ac 1 \
            -Oz > $PREFIX.ZH13.$SAMPLE.SNPs.vcf.gz
    rtg index $PREFIX.ZH13.$SAMPLE.SNPs.vcf.gz
done
```

Extract sample-specific variants from the published set
```
bcftools view -s $SAMPLE --min-ac 1 -Oz Liu.SNPs.chr18.vcf > Liu.SNPs.chr18.$SAMPLE.vcf.gz
rtg index Liu.SNPs.chr18.$SAMPLE.vcf.gz
```

Compare published variants to pggb variants:
```
rtg format \
    -o ZH13.chr18.sdf \
    ZH13.chr18.fa

SAMPLE=C01
for CALLSET in chr18_*/*$SAMPLE*vcf.gz; do
    GRAPH=`echo $CALLSET | sed -E "s/.*chr[0-9][0-9]_(s[0-9]*[0-9])p.*/\1/"`
    PUBSET=../../00_raw/Liu29/variantdata/Liu.SNPs.chr18.$SAMPLE.vcf.gz
    VARTYPE=`echo $CALLSET | sed -E "s/.*\.([^\.]*).vcf.gz/\1/"`
    DIST=100
    # [ Combining cmds with the <() operation does not work for my zsh or whatever other reason ]
    bedtools merge -d $DIST -i $PUBSET >tmp1
    bedtools merge -d $DIST -i $CALLSET >tmp2
    bedtools intersect -a tmp1 -b tmp2 >callable_regions.bed
    rtg vcfeval \
        -b $PUBSET -c $CALLSET \
        -t ZH13.chr18.sdf \
        -T 10 \
        -o vcfeval/$SAMPLE.$GRAPH.$VARTYPE.dist$DIST \
        -e callable_regions.bed
done
```

I tested a few settings, with dist=100/1000, without eval_regions (no dist), and extracting SNPs in Andrea style ('max1') and Joerg style ('SNPs'). Here's the first round output:
```
{echo ""; head -2 <(cat `ls vcfeval/*/summary.txt|head -1`)} | paste - -
for i in vcfeval/*; do echo $i; tail -n1 $i/summary.txt; done | paste - -

        Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure
----------------------------------------------------------------------------------------------------
vcfeval/s30.C01.max1.dist100      None  130223  130223    2410       8814   0.9818   0.9366      0.9587
vcfeval/s30.C01.SNPs              None  130223  130223    41694      35907  0.7575   0.7839      0.7704
vcfeval/s30.C01.SNPs.dist100      None  130223  130223    4281       8886   0.9682   0.9361      0.9519
vcfeval/s30.C01.SNPs.dist1000     None  130223  130223    9020       20949  0.9352   0.8614      0.8968
vcfeval/s50.C01.max1.dist1000     None  129586  129586    9565       21989  0.9313   0.8549      0.8915
vcfeval/s50.C01.max1.dist100      None  129586  129586    3991       9545   0.9701   0.9314      0.9504
vcfeval/s50.C01.SNPs.dist100      None  129586  129586    6076       9656   0.9552   0.9307      0.9428
```

Joerg style keeps more SNPs, and that can also be the reason why Andrea style wins in terms of F-score and alike. Congrats :)

!!! TODO: same for indels
!!! TODO: variant comparison for the Chu 2021 published variants
!!! TODO: for (all) other samples
!!! TODO: for (all) other chromosomes (I would prefer to run only on some chrs because of the runtime, over 24h for one graph w/60 threads)


## Interesting regions

- there is a known duplication associated with resistance to soybean cyst nematode at rhg1 locus on chr18 (https://academic.oup.com/plphys/article/165/2/630/6113099?login=true). Unfortunately, in this set of public lines, it might be that only one line contains a duplication there.
- nice rearrangements leading to gene fusion event on chr19 in E3, a major flowering locus, see Fig. 4 in Liu 2020 paper. How does it look in viz incl. gene annotation?

