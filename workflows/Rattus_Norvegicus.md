# Rattus Norvegicus
We will use the Rattus Norvegicus to build a variation graph and call variants from it.

We use the inbreeding linked-read assemblies of 32 strains from ..., the data was deposited at (...)

Create these folders:
```
mkdir assemblies
mkdir alignments
mkdir parts
```

## Pangenome Sequence Naming
To change the sequence names according to PanSN-spec, we use fastix:

```bash
for f in *_megabubble.fasta.gz; do                                                                                                     
        sample_name=$(echo $f | cut -f 1 -d '_');                                                                                                                    
        fastix -p "${sample_name}#1#" <(zcat $f)> $f._changeid.fa                
        bgzip *_changeid.fa                                                                                                                                             
        samtools faidx *_changeid.fa.gz  
rename "s/.fasta.gz.//" *fa.gz
```

## Sequence divergence
The sequence divergence for each set of chromosomes is high because we are using inbreeding assemblies, so we are using -p 98.

## Sequence partitioning
Partition the assembly contigs by chromosome by mapping each assembly against the reference genome. We use wfmash for the mapping.

```bash
zcat * fasta.gz > rat.fasta.gz
zgrep '_changeid.fasta.gz$' | cut -f 1-2 -d . | sort -V | uniq >haps.list
ref=rn7_supernova_megabubble_changeid.fasta.gz
for hap in $(cat haps.list);
do
    in=assemblies/$hap
    out=alignments/$hap.vs.ref.paf
    wfmash -t 48 -m -N -p 90 -s 20000 $ref $in >$out
done
```
Subset by chromosome:

```bash
(seq 20; echo X; echo Y) | while read i; do awk '$6 ~ "chr'$i'$"' $(ls alignments/*.vs.ref.paf | sort -V) | cut -f 1 | sort -V > parts/chr$i.contigs; done 
(seq 20;echo X;echo Y) | while read i; do xargs samtools faidx rat.fasta.gz < parts/chr$i.contigs > parts/chr$i.pan.fa; done 
(seq 20; echo X; echo Y; echo M) | while read i; do samtools faidx rn7_supernova_megabubble_changeid.fasta.gz rn7#1#chr$i > rn7_changeid.chr$i.fa && cat rn7_changeid.chr$i.fa chr$i.pan.fa > chr$i.pan+ref.fa && bgzip chr$i.pan+ref.fa && samtools faidx chr$i.pan+ref.fa.gz; done
```
This results in chromosome-specific FASTAs in parts/chr*.pan+ref.fa.gz

## Pangenome graph building

We now apply pggb:
```
(seq 20; echo X; echo Y | tr ' ' '\n') | while read i; do sbatch -p workers -c 48 --wrap 'cd /scratch && pggb -i chr'$i'.pan+ref.fa.gz -s 50000 -p 98 -n 32 -k 229 -t 48 -T 32 -O 0.03 -v -Z -V 'rn7:#' -o chr'$i'.pan+ref_rat ; mv /scratch/chr'$i'.pan+ref_rat '$(pwd); done
```

We are considering only the chr1 fr the evaluation.

## Graph statistics
```
odgi stats -i *.smooth.fix.gfa -t 16 -S
#length     nodes     edges       paths
814543152   9523062   13444162    133058
```
## Identify variants with vg
To call variants for each contig, execute:
