# Tomato
Tomato (**Solanum lycoperscium**) is an important fruit crop and a model system for plant biology and breeding.

## Data collection
We use the inbreeding HiFi assemblies from Zhou et.al 2022, the data was deposited at [SolOmics](http://solomics.agis.org.cn/tomato/ftp/genome/)

```
mkdir -p assemblies/tomato
cd assemblies/tomato
cat ../data/tomato.urls| parallel -j 4 'wget -q {} && echo got {}'
```

## Pangenome Sequence Naming
Add prefixs `#1` following the for each assemblies, since all downloaded assemblies are inbreeding lines, most of regions are homozygous.
Remove all the unplaced contigs since high repeat nature of these contigs and can mapped to all 12 chromosomes.

```bash
ls *.fasta.gz | while read f; do
    prefix=${echo $f|cut -f1 -d "."};
    echo $prefix
    fastix -p "${prefix}#1#" < (zcat $f|cut -f 1)|bgzip -@ 32 -c > ${prefix}.fa.gz;
    samtools faidx $prefix.fa.gz
done
```

## Sequence partitioning

All assemblies had scaffolded by `RagTag` by synteny-guided. So we just use all the chromosome information for better graph and parallaize the computing.

```bash
for chr in `seq 1 12`;
do
    for s in `ls *.fa.gz|cut -f1 -d "."`;
    do
        samtools faidx ${s}.fa.gz ${chr} >> chr${chr}.fa
    done
    
    bgzip -@ 24 chr${chr}.fa
    samtools faidx chr${chr}.fa.gz
done
```


## Divergence estimation
Estimated the divergence of each chromosomes to set the parameters of `pggb -p`.

```bash
ls *.fasta.gz | while read f; do mash sketch $f; done
mash triangle *.fa.gz > tomato.mash_triangle.txt
sed 1,1d tomato.mash_triangle.txt | tr '\t' '\n' | grep e -v | sort -g -k 1nr | head -n 5
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
for s in `ls *.fasta.gz|cut -f1 -d "."`;
do 
    bcftools view -s ${s} --threads 4 -v snps Tomato.deepvariant.hifi.vcf.gz|bcftools view -g hom -e 'GT="ref"' -O z -o ${s}.hifi.snps.vcf.gz
    bcftools view -s ${s} --threads 4 -v indels Tomato.deepvariant.hifi.vcf.gz|bcftools view -g hom -e 'GT="ref"' -O z -o ${s}.hifi.indels.vcf.gz
done
```



### Postprocessing the pggb callsets
Filter the variant size and split multiallelic variants

```bash

for v in `ls chr*.fa.gz.*.smooth.fix.SL5.vcf.gz`;
do
    for s in `cat sample.list`;
    do
        bash ../scripts/vcf_preprocess.sh ${v} ${s} 50
    done
done
```

### Benchmarking
`rtg vcfeval` for benchmarking

```bash

for v in `ls chr*.fa.gz.*.smooth.fix.SL5.vcf.gz`;
do
    for s in `cat sample.list`;
    do
        bash ../scripts/small_variant_ ${v} ${s} 50
    done
done

```







