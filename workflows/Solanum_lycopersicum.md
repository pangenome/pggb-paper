# GenusSpecies
Tomato (*Solanum lycoperscium*) is an important fruit crop and a model system for plant biology and breeding.

## Data collection
We use the inbreeding HiFi assemblies from Zhou et.al 2022, the data was deposited at http://solomics.agis.org.cn/tomato/ftp/genome/

```
mkdir -p assemblies/tomato
cd assemblies/tomato
cat ../../docs/data/tomato.urls| parallel -j 4 'wget -q {} && echo got {}'
```

## Pangenome Sequence Naming
Add prefixs `#1` following the for each assemblies, since all downloaded assemblies are inbreeding lines, most of regions are homozygous.
Remove all the unplaced contigs since high repeat nature of these contigs and can mapped to all 12 chromosomes.

```bash
ls *.fasta.gz | while read f; do
    prefix=${echo $f|cut -f1 -d "."};
    echo $prefix
    seqkit fastix -p "${prefix}#1#" < (zcat $f|cut -f 1)|bgzip -@ 32 -c > ${prefix}.fa.gz;
    samtools faidx $prefix.fa.gz
done
```

## Sequence partitioning

Split into different chromosomes for better graph and parallaize the computing.

```bash
for ch in `seq 1 12`;
do
    for i in `ls *.fa.gz|cut -f1 -d "."`;
    do
        samtools faidx ${i}



```


## Divergence estimation
Estimated the divergence of each chromosomes to set the parameters of `pggb -p`.

```bash


```  




## Pangenome graph building

Use the pggb singularity images

```
module load singularity/3.5.2
singularity pull pggb-0406.sif docker://ghcr.io/pangenome/pggb:20220403190741c733e5

for chr in `seq 1 12`;
    do 
        singularity exec pggb-0406.sif pggb -i tomato.${chr}.fa.gz -o tomato.${chr} -t 64 -p 95 -s 10000 -k 49 -n 23 -Z -V "SL5:#"
done

```



## Graph statistics
...


## Small variants evaluation
Call small variants using the DeepVariant using HiFi reads as this repo described [TGG1.1](https://github.com/YaoZhou89/TGG/tree/main/4.Graph_pangenome/1.construction_graph_genome)

### Postprocessing the DeepVariant callsets
Only keep the variant and homozygous variants as baseline of DeppVariant callsets

```bash

```



### Postprocessing the pggb callsets
Filter and 

```bash
```


### Benchmarking





