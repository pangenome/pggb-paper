### All species

The (partial) datasets are in `/lizardfs/guarracino/pggb-paper/sequences`.

Run `pggb` and the Nucmer-based evaluation:

```shell
# Paths
out_folder_graphs=/lizardfs/guarracino/pggb-paper/graphs
out_folder_evaluation=/lizardfs/guarracino/pggb-paper/evaluation
gfa2evaluation=/home/guarracino/pggb-paper/scripts/gfa2evaluation.sh

# Fixed parameters
t=48
POA=asm5
O=0.001


# Homo sapiens
f=/lizardfs/guarracino/pggb-paper/assemblies/hsapiens90.chr6.masked.fa.gz
p=98
s=10000
n=90
k=79
G=700,900,1100
ref=chm13
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.F0001.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c $t -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -t $t -o $out; mv $out $out_folder_graphs"
sbatch -c $t -p 386mem --wrap "hostname; cd /scratch; bash $gfa2evaluation $out_folder_graphs/*.final.gfa $ref $out $t; mv $out $out_folder_evaluation"


# Mus musculus
f=/lizardfs/guarracino/pggb-paper/assemblies/mouse17.chr19.fa.gz
p=98
s=5000
n=17
k=79
G=700,900,1100
ref=GRCm39
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c $t -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -t $t -o $out; mv $out $out_folder_graphs"
sbatch -c $t -p 386mem --wrap "hostname; cd /scratch; bash $gfa2evaluation $out_folder_graphs/*.final.gfa $ref $out $t; mv $out $out_folder_evaluation"


# Rat
f=/lizardfs/guarracino/pggb-paper/assemblies/rat32.chr1.fa.gz
p=98
s=5000
n=32
k=79
G=700,900,1100
ref=rn7
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c $t -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -t $t -o $out; mv $out $out_folder_graphs"
sbatch -c $t -p 386mem --wrap "hostname; cd /scratch; bash $gfa2evaluation $out_folder_graphs/*.final.gfa $ref $out $t; mv $out $out_folder_evaluation"


# Primates
f=/lizardfs/guarracino/pggb-paper/assemblies/primates4.chr6.fa.gz
p=95
s=5000
n=4
k=79
G=700,900,1100
ref=GRC38
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c $t -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -t $t -o $out; mv $out $out_folder_graphs"
sbatch -c $t -p 386mem --wrap "hostname; cd /scratch; bash $gfa2evaluation $out_folder_graphs/*.final.gfa $ref $out $t; mv $out $out_folder_evaluation"

f=/lizardfs/guarracino/pggb-paper/assemblies/primates4.chr7.fa.gz
p=95
s=5000
n=4
k=79
G=700,900,1100
ref=GRC38
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c $t -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -t $t -o $out; mv $out $out_folder_graphs"
sbatch -c $t -p 386mem --wrap "hostname; cd /scratch; bash $gfa2evaluation $out_folder_graphs/*.final.gfa $ref $out $t; mv $out $out_folder_evaluation"


# Arabidopsis thaliana
f=/lizardfs/guarracino/pggb-paper/assemblies/athaliana7.chr1.fa.gz
p=95
s=5000
n=7
k=47
G=700,900,1100
ref=TAIR10
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c $t -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -t $t -o $out; mv $out $out_folder_graphs"
sbatch -c $t -p 386mem --wrap "hostname; cd /scratch; bash $gfa2evaluation $out_folder_graphs/*.final.gfa $ref $out $t; mv $out $out_folder_evaluation"


# Soy
f=/lizardfs/guarracino/pggb-paper/assemblies/soy37.chr18.fa.gz
p=95
s=5000
n=37
k=47
G=700,900,1100
ref=ZH13
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c $t -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -t $t -o $out; mv $out $out_folder_graphs"
sbatch -c $t -p 386mem --wrap "hostname; cd /scratch; bash $gfa2evaluation $out_folder_graphs/*.final.gfa $ref $out $t; mv $out $out_folder_evaluation"


# Tomato
f=/lizardfs/guarracino/pggb-paper/assemblies/tomato23.chr2.fa.gz
p=95
s=5000
n=23
k=47
G=700,900,1100
ref=SL5
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c $t -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -t $t -o $out; mv $out $out_folder_graphs"
sbatch -c $t -p 386mem --wrap "hostname; cd /scratch; bash $gfa2evaluation $out_folder_graphs/*.final.gfa $ref $out $t; mv $out $out_folder_evaluation"


# Yeast
f=/lizardfs/guarracino/pggb-paper/assemblies/scerevisiae8.fa.gz
p=95
s=5000
n=8
k=47
G=700,900,1100
ref=SGDref
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c $t -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -t $t -o $out; mv $out $out_folder_graphs"
sbatch -c $t -p 386mem --wrap "hostname; cd /scratch; bash $gfa2evaluation $out_folder_graphs/*.final.gfa $ref $out $t; mv $out $out_folder_evaluation"


# E. coli
f=/lizardfs/guarracino/pggb-paper/assemblies/ecoli100.fa.gz
p=90
s=5000
n=100
k=47
G=700,900,1100
ref=BH100N_MG2017_3a53c37
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c $t -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -t $t -o $out; mv $out $out_folder_graphs"
sbatch -c $t -p 386mem --wrap "hostname; cd /scratch; bash $gfa2evaluation $out_folder_graphs/*.final.gfa $ref $out $t; mv $out $out_folder_evaluation"
```






TO DELETE:

```shell
echo "Collect statistics"

echo label haplotype region tp.baseline tp.call fp fn precision recall f1.score | tr ' ' '\t' > statistics.tsv

grep '#CHROM' "$PATH_PGGB_VCF" -m 1 | cut -f 10- | tr '\t' '\n' | while read SAMPLE; do
    echo $SAMPLE

    if [[ -f vcfeval/$SAMPLE/easy/summary.txt ]]; then
        # xargs trims whitespaces
        grep None vcfeval/$SAMPLE/easy/summary.txt | tr -s ' ' | xargs | cut -f 2,3,4,5,6,7,8,9 -d ' ' | tr ' ' '\t' | awk -v label=$LABEL -v sample=$SAMPLE -v region="easy" -v OFS='\t' '{print(label, sample, region, $0)}' >> statistics.tsv
        grep None vcfeval/$SAMPLE/hard/summary.txt | tr -s ' ' | xargs | cut -f 2,3,4,5,6,7,8,9 -d ' ' | tr ' ' '\t' | awk -v label=$LABEL -v sample=$SAMPLE -v region="hard" -v OFS='\t' '{print(label, sample, region, $0)}' >> statistics.tsv
    fi
done


# Prepare FPs/FNs for plotting
echo label haplotype region chrom pos ac lv | tr ' ' '\t' > fp.tsv
echo label haplotype region chrom pos ac lv | tr ' ' '\t' > fn.tsv

grep '#CHROM' "$PATH_PGGB_VCF" -m 1 | cut -f 10- | tr '\t' '\n' | while read SAMPLE; do
    echo $SAMPLE

    if [[ -f vcfeval/$SAMPLE/easy/summary.txt ]]; then
        PATH_FP_VCF=vcfeval/$SAMPLE/easy/fp.vcf.gz
        zcat $PATH_FP_VCF | cut -f -8 | vcf2tsv | cut -f 1,2,8,13 | awk -v label=$LABEL -v sample=$SAMPLE -v region="easy" -v OFS='\t' 'NR > 1 { print label, sample, region, $0 }' >> fp.tsv

        PATH_FP_VCF=vcfeval/$SAMPLE/hard/fp.vcf.gz
        zcat $PATH_FP_VCF | cut -f -8 | vcf2tsv | cut -f 1,2,8,13 | awk -v label=$LABEL -v sample=$SAMPLE -v region="hard" -v OFS='\t' 'NR > 1 { print label, sample, region, $0 }' >> fp.tsv


        PATH_FN_VCF=vcfeval/$SAMPLE/easy/fn.vcf.gz
        zcat $PATH_FN_VCF | cut -f -8 | vcf2tsv | cut -f 1,2,8,13 | awk -v label=$LABEL -v sample=$SAMPLE -v region="easy" -v OFS='\t' 'NR > 1 { print label, sample, region, $0 }' >> fn.tsv

        PATH_FN_VCF=vcfeval/$SAMPLE/hard/fn.vcf.gz
        zcat $PATH_FN_VCF | cut -f -8 | vcf2tsv | cut -f 1,2,8,13 | awk -v label=$LABEL -v sample=$SAMPLE -v region="hard" -v OFS='\t' 'NR > 1 { print label, sample, region, $0 }' >> fn.tsv
    fi
done
```
