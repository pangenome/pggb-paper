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
RUN_WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-ad8aebae1be96847839778af534866bc9545adb9
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


```shell
mkdir -p /lizardfs/guarracino/pggb-paper/hsapiens/graphs

RUN_PGGB=/home/guarracino/tools/pggb/pggb-a4a6668d9ece42c80ce69dc354f0cb59a849286f

( echo 1 16 2 3 4 5 6 7 8 X 9 10 11 12 13 14 15 17 18 19 20 21 22 Y | tr ' ' '\n') | while read i; do
  sbatch -p workers -c 48 --wrap "hostname; cd /scratch && $RUN_PGGB -i /lizardfs/guarracino/pggb-paper/hsapiens/assemblies/parts/chr$i.pan.fa.gz -o chr$i.pan -t 48 -p 98 -s 50000 -n 8 -k 311 -O 0.03 -T 48 ; mv /scratch/chr$i.pan /lizardfs/guarracino/pggb-paper/hsapiens/graphs/"
done
```


## Graph statistics
...


## Small variants evaluation
...

