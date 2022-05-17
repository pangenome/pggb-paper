# HomoSapiens


## Data collection

Get the URLs of the assemblies:

```shell
mkdir -p /lizardfs/guarracino/pggb-paper/hsapiens/assemblies
cd /lizardfs/guarracino/pggb-paper/hsapiens/assemblies

wget https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/assembly_index/Year1_assemblies_v2_genbank.index
grep 'chm13\|h38' Year1_assemblies_v2_genbank.index | awk '{ print $2 }' | sed 's%s3://human-pangenomics/working/%https://s3-us-west-2.amazonaws.com/human-pangenomics/working/%g' > refs.urls
grep 'chm13\|h38' -v Year1_assemblies_v2_genbank.index  | grep HPRC_PLUS | awk '{ print $2; print $3 }' | sed 's%s3://human-pangenomics/working/%https://s3-us-west-2.amazonaws.com/human-pangenomics/working/%g' > samples.urls
```

Download 6 haplotypes (3 diploid samples) and the references:

```shell
cat refs.urls <(grep HG03 samples.urls | head -n 6) | parallel -j 4 'wget -q {} && echo got {}'
```

Unpack the assemblies:

```shell
ls *.f1_assembly_v2_genbank.fa.gz | while read f; do echo $f; gunzip $f && samtools faidx $(basename $f .gz); done
```

## Pangenome Sequence Naming

Add a prefix to the reference sequences:

```shell
fastix  -p 'grch38#' <(zcat GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz) | bgzip -c -@ 48 > grch38_full.fa.gz && samtools faidx grch38_full.fa.gz
fastix  -p 'chm13#' <(zcat chm13.draft_v1.1.fasta.gz) | bgzip -c -@ 48 > chm13.fa.gz && samtools faidx chm13.fa.gz

# Remove unplaced contigs from grch38 that are (hopefully) represented in chm13
samtools faidx grch38_full.fa.gz $(cat grch38_full.fa.gz.fai | cut -f 1 | grep -v _ ) | bgzip -@ 48 -c > grch38.fa.gz

# Cleaning
rm chm13.draft_v1.1.fasta.gz GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz grch38_full.fa.gz.*
```

## Sequence partitioning

Prepare the references:

```shell
zcat chm13.fa.gz grch38.fa.gz | bgzip -c -@ 48 > chm13+grch38.fa.gz && samtools faidx chm13+grch38.fa.gz
```

We partition contigs by chromosome by mapping each assembly against the scaffolded references:

```shell
RUN_WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-48493fbe2d7e65b7e371b8940590c88dac13e4ad
PATH_REFERENCE_FA_GZ=/lizardfs/guarracino/pggb-paper/hsapiens/assemblies/chm13+grch38.fa.gz

mkdir -p /lizardfs/guarracino/pggb-paper/hsapiens/assemblies/partitioning

ls *.f1_assembly_v2_genbank.fa | while read FASTA; do
  NAME=$(basename $FASTA .fa);
  echo $NAME
  
  PATH_PAF=/lizardfs/guarracino/pggb-paper/hsapiens/assemblies/partitioning/$NAME.vs.ref.paf
  sbatch -p headnode -c 12 --wrap "$RUN_WFMASH $PATH_REFERENCE_FA_GZ /lizardfs/guarracino/pggb-paper/hsapiens/assemblies/$FASTA -s 50k -p 90 -N -m -t 12 > $PATH_PAF"
done
```

Collect unmapped contigs and remap them in split mode:

```shell
ls *.f1_assembly_v2_genbank.fa | while read FASTA; do
  NAME=$(basename $FASTA .fa);
  echo $NAME
  
  PATH_PAF=/lizardfs/guarracino/pggb-paper/hsapiens/assemblies/partitioning/$NAME.vs.ref.paf
  comm -23 <(cut -f 1 $FASTA.fai | sort) <(cut -f 1 $PATH_PAF | sort) > $NAME.unaligned.txt
  
  wc -l $NAME.unaligned.txt
  if [[ $(wc -l $NAME.unaligned.txt | cut -f 1 -d\ ) != 0 ]];
  then 
    samtools faidx $FASTA $(tr '\n' ' ' < $NAME.unaligned.txt) > $NAME.unaligned.fa
    samtools faidx $NAME.unaligned.fa
    
    PATH_NO_SPLIT_PAF=/lizardfs/guarracino/pggb-paper/hsapiens/assemblies/partitioning/$NAME.vs.ref.no_split.paf
    sbatch -p headnode -c 12 --wrap "wfmash $PATH_REFERENCE_FA_GZ /lizardfs/guarracino/pggb-paper/hsapiens/assemblies/$NAME.unaligned.fa -s 50k -p 90 -m -t 12 > $PATH_NO_SPLIT_PAF"
  fi
done
```

Collect our best mapping for each of our attempted split rescues:

```shell
cd /lizardfs/guarracino/pggb-paper/hsapiens/assemblies/partitioning

ls *.vs.ref.no_split.paf | while read PAF; do
  cat $PAF | awk -v OFS='\t' '{ print $1,$11,$0 }' | sort -n -r -k 1,2 | \
    awk -v OFS='\t' '$1 != last { print($3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15); last = $1; }'
done > rescues.paf
```

Subset by chromosome, including the references:

```shell
cd /lizardfs/guarracino/pggb-paper/hsapiens/assemblies

DIR_PARTITIONING=/lizardfs/guarracino/pggb-paper/hsapiens/assemblies/partitioning

mkdir -p parts

( seq 22; echo X; echo Y; echo M ) | while read i; do
  awk '$6 ~ "chr'$i'$"' $(ls $DIR_PARTITIONING/*.vs.ref.paf | \
    grep -v unaligned | sort; echo $DIR_PARTITIONING/rescues.paf) | cut -f 1 | sort | uniq \
    > parts/chr$i.contigs;
done

( seq 22; echo X; echo Y; echo M ) | while read i; do
  echo chr$i
  samtools faidx chm13+grch38.fa.gz chm13#chr$i grch38#chr$i > parts/chr$i.pan.fa
  
  ls *.f1_assembly_v2_genbank.fa | while read FASTA; do
    NAME=$(basename $FASTA .fa);
    echo chr$i $NAME
    
    samtools faidx $FASTA $( comm -12 <(cut -f 1 $FASTA.fai | sort) <(sort parts/chr$i.contigs) ) >> parts/chr$i.pan.fa
  done
  
  bgzip -@ 48 parts/chr$i.pan.fa && samtools faidx parts/chr$i.pan.fa.gz
done
```

## Divergence estimation

The assemblies are at the contig level, so we can't just estimate the divergence with `mash`.
Since human presents a low sequence divergence, so we will set the mapping identity (`-p` parameter`) in `pggb` to `98`.


## Pangenome graph building

Build the pangenome graph for each chromosome:

```shell
mkdir -p /lizardfs/guarracino/pggb-paper/hsapiens/graphs

RUN_PGGB=/home/guarracino/tools/pggb/pggb-4f71d69f4c0bcacbd9b286c1aaed86e7c37fc0a2

( echo 1 16 2 3 4 5 6 7 8 X 9 10 11 12 13 14 15 17 18 19 20 21 22 Y | tr ' ' '\n') | while read i; do
  sbatch -p workers -c 48 --wrap "hostname; cd /scratch && $RUN_PGGB -i /lizardfs/guarracino/pggb-paper/hsapiens/assemblies/parts/chr$i.pan.fa.gz -o chr$i.pan -t 48 -p 98 -s 50000 -n 8 -k 311 -O 0.03 -T 48 ; mv /scratch/chr$i.pan /lizardfs/guarracino/pggb-paper/hsapiens/graphs/"
done
```

## Graph statistics

To collect basic graph statistics, execute:

```shell
RUN_ODGI=/home/guarracino/tools/odgi/bin/odgi-4d4acae6d7e27fed7cdbfc5c23c3e32d28092e89

( echo 1 16 2 3 4 5 6 7 8 X 9 10 11 12 13 14 15 17 18 19 20 21 22 Y | tr ' ' '\n') | while read i; do
    sbatch -p workers -c 48 --wrap "hostname; $RUN_ODGI stats -i /lizardfs/guarracino/pggb-paper/hsapiens/graphs/chr$i.pan/chr$i.pan.fa.gz.a6b8705.04f1c29.65f7a81.smooth.final.og -S > /lizardfs/guarracino/pggb-paper/hsapiens/graphs/chr$i.pan/chr$i.pan.fa.gz.a6b8705.04f1c29.65f7a81.smooth.final.stats.tsv"
done
```


## Small variants evaluation

### Identify variants with vg

Call variants for each haplotype:

```shell
# TO FIX: vg deconstruct does not support that mode yet
# Call variants for each haplotype
( echo 1 16 2 3 4 5 6 7 8 X 9 10 11 12 13 14 15 17 18 19 20 21 22 Y | tr ' ' '\n') | while read i; do
  PATH_GFA=/lizardfs/guarracino/pggb-paper/hsapiens/graphs/chr$i.pan/chr$i.pan.fa.gz.a6b8705.04f1c29.65f7a81.smooth.final.gfa
  PATH_SED_GFA=/lizardfs/guarracino/pggb-paper/hsapiens/graphs/chr$i.pan/chr$i.pan.fa.gz.a6b8705.04f1c29.65f7a81.smooth.final.sed.gfa
  # sed replaces only the first instance on a line by default (without the /g modifier)
  # To have names like NA21309-1#1#JAHEPC010000450.1 and call haploid genotypes with -H
  sed 's/#/-/' $PATH_GFA | sed 's/#/#1#/' > $PATH_SED_GFA

  sbatch -p workers -c 48 --wrap "hostname; cd /scratch && \time -v vg deconstruct -P chm13 -H '#' -e -a -t 48 $PATH_SED_GFA > /lizardfs/guarracino/pggb-paper/hsapiens/graphs/chr$i.pan/chr$i.pan.fa.gz.a6b8705.04f1c29.65f7a81.smooth.final.chm13.haplo.vcf"
done

# Revert names in the VCF files
( echo 1 16 2 3 4 5 6 7 8 X 9 10 11 12 13 14 15 17 18 19 20 21 22 Y | tr ' ' '\n') | while read i; do
  echo chr$i
  
  PATH_VCF=/lizardfs/guarracino/pggb-paper/hsapiens/graphs/chr$i.pan/chr$i.pan.fa.gz.a6b8705.04f1c29.65f7a81.smooth.final.chm13.haplo.vcf
  grep '^##' $PATH_VCF | sed 's/chm13-/chm13#/g' > x.vcf
  grep '^#CHROM' $PATH_VCF | sed 's/-/#/g' >> x.vcf
  grep '^#' $PATH_VCF -v | sed 's/chm13-/chm13#/g' >> x.vcf
  mv x.vcf $PATH_VCF
done
```

Take SNPs from the PGGB VCF files:

```shell
( echo 1 16 2 3 4 5 6 7 8 X 9 10 11 12 13 14 15 17 18 19 20 21 22 | tr ' ' '\n') | while read i; do
  echo chr$i
  
  PATH_CHR_REF_FA=/lizardfs/guarracino/pggb-paper/hsapiens/nucmer/chm13#chr$i.fa  
  PATH_VCF=/lizardfs/guarracino/pggb-paper/hsapiens/graphs/chr$i.pan/chr$i.pan.fa.gz.a6b8705.04f1c29.65f7a81.smooth.final.chm13.haplo.vcf
  grep '#CHROM' $PATH_VCF -m 1 | cut -f 10- | tr '\t' '\n' | while read HAPLO; do
    echo $HAPLO
    bash /home/guarracino/tools/pggb/scripts/vcf_preprocess.sh \
        $PATH_VCF \
        $HAPLO \
        1 \
        $PATH_CHR_REF_FA
  done
done
```

### Identify variants with `nucmer`

Align each contig against the reference with `nucmer`:

```shell
mkdir - p /lizardfs/guarracino/pggb-paper/hsapiens/nucmer
cd /lizardfs/guarracino/pggb-paper/hsapiens/nucmer

( echo 1 16 2 3 4 5 6 7 8 X 9 10 11 12 13 14 15 17 18 19 20 21 22 | tr ' ' '\n') | while read i; do
  echo chr$i
  
  PATH_CHR_FA_GZ=/lizardfs/guarracino/pggb-paper/hsapiens/assemblies/parts/chr$i.pan.fa.gz
  PATH_CHR_REF_FA=/lizardfs/guarracino/pggb-paper/hsapiens/nucmer/chm13#chr$i.fa
  
  samtools faidx $PATH_CHR_FA_GZ "chm13#chr$i" > $PATH_CHR_REF_FA

  cut -f 1 $PATH_CHR_FA_GZ.fai | grep chm13 -v | while read CONTIG; do   
    PREFIX=/lizardfs/guarracino/pggb-paper/hsapiens/nucmer/$CONTIG
    sbatch -p workers -c 6 --wrap "hostname; samtools faidx $PATH_CHR_FA_GZ $CONTIG > $PREFIX.fa && \time -v nucmer $PATH_CHR_REF_FA $PREFIX.fa --prefix $PREFIX -t 6"
  done
done
```

Using the `nucmer2vcf.R` script, generate VCF files for each sequence with respect to the reference:

```shell
NUCMER_VERSION="4.0.0beta2"

( echo 1 16 2 3 4 5 6 7 8 X 9 10 11 12 13 14 15 17 18 19 20 21 22 | tr ' ' '\n') | while read i; do
  PATH_CHR_FA_GZ=/lizardfs/guarracino/pggb-paper/hsapiens/assemblies/parts/chr$i.pan.fa.gz
  PATH_CHR_REF_FA=/lizardfs/guarracino/pggb-paper/hsapiens/nucmer/chm13#chr$i.fa
  
  cut -f 1 $PATH_CHR_FA_GZ.fai | grep chm13 -v | while read CONTIG; do
    echo $CONTIG
    
    PREFIX=/lizardfs/guarracino/pggb-paper/hsapiens/nucmer/$CONTIG
    show-snps -THC $PREFIX.delta > $PREFIX.var.txt
    #show-snps -TH $PREFIX.delta | cut -f 1-6,9-12 > $PREFIX.var.txt # For taking also variants from ambiguous alignments
    
    # Check if there are variants
    if [[ $(wc -l $PREFIX.var.txt | cut -f 1 -d\ ) != 0 ]]; then
      Rscript /home/guarracino/tools/pggb/scripts/nucmer2vcf.R $PREFIX.var.txt $CONTIG $PATH_CHR_REF_FA $NUCMER_VERSION $PREFIX.vcf
      bgzip -@ 48 $PREFIX.vcf
      tabix $PREFIX.vcf.gz
    fi
  done
done
```

Merge variants by haplotype:

```shell
cd /lizardfs/guarracino/pggb-paper/hsapiens/nucmer

REF_PREFIX=chm

( echo 1 16 2 3 4 5 6 7 8 X 9 10 11 12 13 14 15 17 18 19 20 21 22 | tr ' ' '\n') | while read i; do
  mkdir -p /lizardfs/guarracino/pggb-paper/hsapiens/nucmer/haplotype/chr$i
  
  PATH_CHR_FA_GZ=/lizardfs/guarracino/pggb-paper/hsapiens/assemblies/parts/chr$i.pan.fa.gz

  cut -f 1 $PATH_CHR_FA_GZ.fai | grep chm13 -v | cut -f 1,2 -d '#' | sort | uniq | while read HAPLO; do
    echo $HAPLO
    
    grep "^$HAPLO" $PATH_CHR_FA_GZ.fai | cut -f 1 | while read CONTIG; do
      PATH_VCF=$CONTIG.vcf.gz
      if [[ -f $PATH_VCF ]]; then
        # Check if there are valid variants
        if [[ $(zgrep '^#' $PATH_VCF -v | grep "^$REF_PREFIX" | wc -l) != 0 ]]; then
          echo $PATH_VCF >> x.txt
        fi
      fi      
    done
    
    bcftools concat \
      --file-list x.txt \
      --naive-force | \
      bcftools reheader --samples <(echo $HAPLO) > haplotype/chr$i/$HAPLO.unsorted.vcf.gz
    rm x.txt
      
    #bcftools sort $HAPLO.unsorted.vcf.gz > $HAPLO.vcf.gz # It doesn't work
    zgrep "^#" haplotype/chr$i/$HAPLO.unsorted.vcf.gz > haplotype/chr$i/$HAPLO.vcf
    zgrep -v "^#" haplotype/chr$i/$HAPLO.unsorted.vcf.gz | sort -k 1,1V -k 2,2n >> haplotype/chr$i/$HAPLO.vcf
    bgzip -@ 48 haplotype/chr$i/$HAPLO.vcf
    tabix haplotype/chr$i/$HAPLO.vcf.gz
      
    rm haplotype/chr$i/$HAPLO.unsorted.vcf.gz
  done
done
```

### Evaluation

Prepare the reference in `SDF` format for variant evaluation with `rtg vcfeval`:

```shell
RUN_RTG=/gnu/store/3vmp4dw8y0r49h0hbjbgv3bckgvz4k0m-rtg-tools-3.11/rtg

( echo 1 16 2 3 4 5 6 7 8 X 9 10 11 12 13 14 15 17 18 19 20 21 22 | tr ' ' '\n') | while read i; do
  PATH_CHR_REF_FA=/lizardfs/guarracino/pggb-paper/hsapiens/nucmer/chm13#chr$i.fa
  
  $RUN_RTG format -o /lizardfs/guarracino/pggb-paper/hsapiens/nucmer/chm13#chr$i.sdf /lizardfs/guarracino/pggb-paper/hsapiens/nucmer/chm13#chr$i.fa
done
```

Compare nucmer-based SNPs with PGGB-based SNPs:

```shell
( echo 1 16 2 3 4 5 6 7 8 X 9 10 11 12 13 14 15 17 18 19 20 21 22 | tr ' ' '\n') | while read i; do
  PATH_CHR_REF_SDF=/lizardfs/guarracino/pggb-paper/hsapiens/nucmer/chm13#chr$i.sdf
  
  PATH_CHR_FA_GZ=/lizardfs/guarracino/pggb-paper/hsapiens/assemblies/parts/chr$i.pan.fa.gz
  cut -f 1 $PATH_CHR_FA_GZ.fai | grep chm13 -v | cut -f 1,2 -d '#' | sort | uniq | while read HAPLO; do
    echo $HAPLO
    
    PATH_NUCMER_VCF=/lizardfs/guarracino/pggb-paper/hsapiens/nucmer/haplotype/chr$i/$HAPLO.vcf.gz
    PATH_PGGB_VCF=/lizardfs/guarracino/pggb-paper/hsapiens/graphs/chr$i.pan/chr$i.pan.fa.gz.a6b8705.04f1c29.65f7a81.smooth.final.chm13.haplo.vcf.$HAPLO.max1.vcf.gz

    # Merge regions closer than 1000 bps to define the callable regions where to evaluate the variants
    dist=1000

    $RUN_RTG vcfeval \
        -t $PATH_CHR_REF_SDF \
        -b $PATH_NUCMER_VCF \
        -c $PATH_PGGB_VCF \
        -T 16 \
        -e <(bedtools intersect -a <(bedtools merge -d $dist -i $PATH_NUCMER_VCF ) -b <(bedtools merge -d $dist -i $PATH_PGGB_VCF)) \
        -o vcfeval/chr$i/$HAPLO
  done
done
```

Collect statistics:

```shell
echo chromosome haplotype precision recall f1.score | tr ' ' '\t' > /lizardfs/guarracino/pggb-paper/hsapiens/vcfeval/statistics.tsv
( echo 1 16 2 3 4 5 6 7 8 X 9 10 11 12 13 14 15 17 18 19 20 21 22 | tr ' ' '\n') | while read i; do
    cd /lizardfs/guarracino/pggb-paper/hsapiens/vcfeval/chr$i
    
    # TO FIX: sed 's/grch38#chr[0-9X]*/grch38#1/g' is necessary because the references do not respect the PanSN-spec convention
    grep None */summary.txt | sed 's,/summary.txt:,,' | tr -s ' ' | cut -f 1,7,8,9 -d ' ' | tr ' ' '\t' | sed 's/grch38#chr[0-9X]*/grch38#1/g' | awk -v chr=chr$i -v OFS='\t' '{print(chr, $0)}' >> /lizardfs/guarracino/pggb-paper/hsapiens/vcfeval/statistics.tsv
done
```
