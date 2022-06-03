# Mus musculus
We will use the Mus musculus to build a variation graph and call variants from it.

We use the inbreeding linked-read assemblies of 148 extant BXD strains from D.G. Ashbrook et.al 2022, the data was deposited at (...)

Create these folders:
```
mkdir assemblies
mkdir alignments
mkdir parts
```

## Pangenome Sequence Naming
To change the sequence names according to PanSN-spec, we use fastix:

```bash
for f in *.fasta.gz; do
        sample_name=$(echo $f | cut -f 1 -d '_');
        fastix -p "${sample_name}#1#" <(zcat $f)> /home/flaviav/data/BXD_pggb_paper/assemblies/$f._changeid.fa
        bgzip *_changeid.fa
        samtools faidx changeid.fa.gz
rename "s/.fasta.gz.//" *fa.gz
```

## Sequence divergence
The sequence divergence for each set of chromosomes is hight because we are using inbreeding assemblies, so we are using -p 98.

## Sequence partitioning
Partition the assembly contigs by chromosome by mapping each assembly against the reference genome. We use wfmash for the mapping.

```bash
zcat *fasta.gz > mouse.fasta.gz
ls *fasta.gz| sort -V | uniq >haps.list
ref=UCSC_mm10.changeid.fa.gz
for hap in $(cat haps.list);
do
    in=assemblies/$hap
    out=alignments/$hap.vs.ref.paf
    wfmash -t 48 -m -N -p 90 -s 20000 $ref $in >$out
done
```
Subset by chromosome:

```bash
(seq 19; echo X; echo Y) | while read i; do awk '$6 ~ "chr'$i'$"' $(ls alignments/*.vs.ref.paf | sort -V) | cut -f 1 | sort -V > parts/chr$i.contigs; done 
(seq 19;echo X;echo Y) | while read i; do xargs samtools faidx mouse.fasta.gz < parts/chr$i.contigs > parts/chr$i.pan.fa; done
(seq 19; echo X; echo Y) | while read i; do samtools faidx UCSC_mm10.changeid.fa REF#1#chr$i > UCSC_mm10.changeid.chr$i.fa && cat UCSC_mm10.changeid.chr$i.fa chr$i.pan.fa > chr$i.pan+ref.fa && bgzip chr$i.pan+ref.fa && samtools faidx chr$i.pan+ref.fa.gz; done
```

This results in chromosome-specific FASTAs in parts/chr*.pan+ref.fa.gz

## Pangenome graph building

We now apply pggb:
```
( seq 19; echo X; echo Y | tr ' ' '\n') | while read i; do sbatch -p workers -c 48 --wrap 'cd /scratch && /gnu/store/ag0rfhmfjdbmj612mn0xc7kiacf12j8v-pggb-0.3.0+f80557a-2/bin/pggb -i /lizardfs/flaviav/mouse/pggb_paper/parts/chr'$i'.pan+ref.fa.gz -s 50000 -p 98 -n 148 -k 229 -t 48 -T 32 -O 0.03 -v -Z -V 'REF:#' -o chr'$i'.pan+ref_mouse ; mv /scratch/chr'$i'.pan+ref_mouse '$(pwd); done
```

We are considering only the chr19 (61,420,004 bp) for the evaluation.


## Graph statistics
```
odgi stats -i *.smooth.fix.gfa -t 16 -S
#length     nodes   edges   paths
355179727   5034423 8112765 73092
```
## Identify variants with vg

To call variants for each contig, execute:
