### All species

The (partial) datasets are in `/lizardfs/guarracino/pggb-paper/sequences`.

`pggb` command lines:

```shell
out_folder=/lizardfs/guarracino/pggb-paper/graphs

# Fixed parameters
t=48
POA=asm5
O=0.03

# Homo sapiens
f=/lizardfs/guarracino/pggb-paper/sequences/hsapiens90.chr6.masked.fa.gz
p=98
s=10000
n=90
F=0.001
k=79
G=4001,4507
ref=chm13
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.F0001.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c 48 -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -F $F -k $k -P $POA -O $O -G $G -V $ref:# -t $t -o $out; mv $out $out_folder"

f=/lizardfs/guarracino/pggb-paper/sequences/hsapiens90.chr6.masked.fa.gz
p=98
s=5000
n=90
k=79
G=4001,4507
ref=chm13
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c 48 -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -V $ref:# -t $t -o $out; mv $out $out_folder"

f=/lizardfs/guarracino/pggb-paper/sequences/hsapiens90.chr6.masked.fa.gz
p=98
s=5000
n=90
k=271
G=4001,4507
ref=chm13
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c 48 -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -V $ref:# -t $t -o $out; mv $out $out_folder"

# Mus musculus
f=/lizardfs/guarracino/pggb-paper/sequences/mouse17.chr19.fa.gz
p=98
s=5000
n=17
k=79
G=4001,4507
ref=GRCm39
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c 48 -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -V $ref:# -t $t -o $out; mv $out $out_folder"

# Rat
f=/lizardfs/guarracino/pggb-paper/sequences/rat32.chr1.fa.gz
p=98
s=5000
n=32
k=271
G=4001,4507
ref=rn7
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c 48 -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -V $ref:# -t $t -o $out; mv $out $out_folder"

# Primates
f=/lizardfs/guarracino/pggb-paper/sequences/primates4.chr6.fa.gz
p=95
s=5000
n=4
k=47
G=4001,4507
ref=GRC38
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c 48 -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -V $ref:# -t $t -o $out; mv $out $out_folder"

f=/lizardfs/guarracino/pggb-paper/sequences/primates4.chr7.fa.gz
p=95
s=5000
n=4
k=47
G=4001,4507
ref=GRC38
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c 48 -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -V $ref:# -t $t -o $out; mv $out $out_folder"

# Arabidopsis thaliana
f=/lizardfs/guarracino/pggb-paper/sequences/athaliana7.chr1.fa.gz
p=95
s=5000
n=7
k=47
G=4001,4507
ref=TAIR10
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c 48 -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -V $ref:# -t $t -o $out; mv $out $out_folder"

# Soy
f=/lizardfs/guarracino/pggb-paper/sequences/soy37.chr18.fa.gz
p=95
s=5000
n=37
k=47
G=4001,4507
ref=ZH13
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c 48 -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -V $ref:# -t $t -o $out; mv $out $out_folder"

# Tomato
f=/lizardfs/guarracino/pggb-paper/sequences/tomato23.chr2.fa.gz
p=95
s=5000
n=23
k=47
G=4001,4507
ref=SL5
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c 48 -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -V $ref:# -t $t -o $out; mv $out $out_folder"

# Yeast
f=/lizardfs/guarracino/pggb-paper/sequences/scerevisiae8.fa.gz
p=95
s=5000
n=8
k=47
G=4001,4507
ref=SGDref
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c 48 -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -V $ref:# -t $t -o $out; mv $out $out_folder"

# E. coli
f=/lizardfs/guarracino/pggb-paper/sequences/ecoli100.fa.gz
p=90
s=5000
n=100
k=47
G=4001,4507
ref=BH100N_MG2017_3a53c37
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c 48 -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -V $ref:# -t $t -o $out; mv $out $out_folder"
```

t=48
POA=asm5
O=0.03
f=/lizardfs/flaviav/rat/pggb_paper/parts/chr12.pan+ref.fa.gz
p=98
s=1000
n=32
k=79
G=4001,4507
ref=rn7
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c 48 -p workers --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -V $ref:# -t $t -o $out; mv $out $out_folder"

t=48
POA=asm5
O=0.03
f=/lizardfs/flaviav/rat/pggb_paper/parts/chr12.pan+ref.fa.gz
p=98
s=1000
n=32
k=79
G=4001,4507
ref=rn7
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.F0001.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c 48 -p workers --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -F 0.001 -k $k -P $POA -O $O -G $G -V $ref:# -t $t -o $out; mv $out $out_folder"

t=48
POA=asm5
O=0.03
f=/lizardfs/flaviav/rat/pggb_paper/parts/chr12.pan+ref.fa.gz
p=98
s=5000
n=32
k=79
G=4001,4507
ref=rn7
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c 48 -p workers --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -V $ref:# -t $t -o $out; mv $out $out_folder"

Nucmer-based evaluation:


```shell
PATH_VCF_PREPROCESS=/home/guarracino/tools/pggb/scripts/vcf_preprocess.sh
PATH_NUCMER_2_VCF=/home/guarracino/tools/pggb/scripts/nucmer2vcf.R


sbatch -c 48 -p workers --wrap "hostname; cd /scratch && bash /home/guarracino/pggb-paper/scripts/gfa2evaluation.sh /lizardfs/guarracino/pggb-paper/graphs/hsapiens90.chr6_p98.s1000.n90.k271.G13117-13219.chm13/*.final.gfa chm13 /lizardfs/guarracino/pggb-paper/evaluation/hsapiens90.chr6_p98.s1000.n90.k271.G13117-13219.chm13/ $PATH_VCF_PREPROCESS $PATH_NUCMER_2_VCF 48"

sbatch -c 48 -p workers --wrap "hostname; cd /scratch && bash /home/guarracino/pggb-paper/scripts/gfa2evaluation.sh /lizardfs/guarracino/pggb-paper/graphs/hsapiens90.chr6_p98.s5000.n90.k271.G13117-13219.chm13/*.final.gfa chm13 /lizardfs/guarracino/pggb-paper/evaluation/hsapiens90.chr6_p98.s5000.n90.k271.G13117-13219.chm13/ $PATH_VCF_PREPROCESS $PATH_NUCMER_2_VCF 48"

sbatch -c 48 -p workers --wrap "hostname; cd /scratch && bash /home/guarracino/pggb-paper/scripts/gfa2evaluation.sh /lizardfs/guarracino/pggb-paper/graphs/hsapiens90.chr6_p98.s5000.n90.k271.G4001-4507.chm13/*.final.gfa chm13 /lizardfs/guarracino/pggb-paper/evaluation/hsapiens90.chr6_p98.s5000.n90.k271.G4001-4507.chm13/ $PATH_VCF_PREPROCESS $PATH_NUCMER_2_VCF 48"


sbatch -c 48 -p workers --wrap "hostname; cd /scratch && bash /home/guarracino/pggb-paper/scripts/gfa2evaluation.sh /lizardfs/guarracino/pggb-paper/graphs/mouse17.chr19_p98.s5000.n17.k79.G4001-4507.GRCm39/*.final.gfa GRCm39 /lizardfs/guarracino/pggb-paper/evaluation/mouse17.chr19_p98.s5000.n17.k79.G4001-4507.GRCm39/ $PATH_VCF_PREPROCESS $PATH_NUCMER_2_VCF 48"


sbatch -c 48 -p workers --wrap "hostname; cd /scratch && bash /home/guarracino/pggb-paper/scripts/gfa2evaluation.sh /lizardfs/guarracino/pggb-paper/graphs/rat32.chr1_p98.s1000.n32.k271.G4001-4507.rn7/*.final.gfa rn7 /lizardfs/guarracino/pggb-paper/evaluation/rat32.chr1_p98.s1000.n32.k271.G4001-4507.rn7/ $PATH_VCF_PREPROCESS $PATH_NUCMER_2_VCF 48"

sbatch -c 48 -p workers --wrap "hostname; cd /scratch && bash /home/guarracino/pggb-paper/scripts/gfa2evaluation.sh /lizardfs/guarracino/pggb-paper/graphs/primates4.chr6_p95.s5000.n4.k47.G4001-4507.GRC38/*.final.gfa GRC38 /lizardfs/guarracino/pggb-paper/evaluation/primates4.chr6_p95.s5000.n4.k47.G4001-4507.GRC38/ $PATH_VCF_PREPROCESS $PATH_NUCMER_2_VCF 48"

sbatch -c 48 -p workers --wrap "hostname; cd /scratch && bash /home/guarracino/pggb-paper/scripts/gfa2evaluation.sh /lizardfs/guarracino/pggb-paper/graphs/primates4.chr7_p95.s5000.n4.k47.G4001-4507.GRC38/*.final.gfa GRC38 /lizardfs/guarracino/pggb-paper/evaluation/primates4.chr7_p95.s5000.n4.k47.G4001-4507.GRC38/ $PATH_VCF_PREPROCESS $PATH_NUCMER_2_VCF 48"

sbatch -c 48 -p workers --wrap "hostname; cd /scratch && bash /home/guarracino/pggb-paper/scripts/gfa2evaluation.sh /lizardfs/guarracino/pggb-paper/graphs/athaliana7.chr1_p95.s5000.n7.k47.G4001-4507.TAIR10/*.final.gfa TAIR10 /lizardfs/guarracino/pggb-paper/evaluation/athaliana7.chr1_p95.s5000.n7.k47.G4001-4507.TAIR10/ $PATH_VCF_PREPROCESS $PATH_NUCMER_2_VCF 48"

sbatch -c 48 -p workers --wrap "hostname; cd /scratch && bash /home/guarracino/pggb-paper/scripts/gfa2evaluation.sh /lizardfs/guarracino/pggb-paper/graphs/soy37.chr18_p95.s5000.n37.k47.G4001-4507.ZH13/*.final.gfa ZH13 /lizardfs/guarracino/pggb-paper/evaluation/soy37.chr18_p95.s5000.n37.k47.G4001-4507.ZH13 $PATH_VCF_PREPROCESS $PATH_NUCMER_2_VCF 48"

sbatch -c 48 -p workers --wrap "hostname; cd /scratch && bash /home/guarracino/pggb-paper/scripts/gfa2evaluation.sh /lizardfs/guarracino/pggb-paper/graphs/tomato23.chr2_p95.s5000.n23.k47.G4001-4507.SL5/*.final.gfa SL5 /lizardfs/guarracino/pggb-paper/evaluation/tomato23.chr2_p95.s5000.n23.k47.G4001-4507.SL5/ $PATH_VCF_PREPROCESS $PATH_NUCMER_2_VCF 48"


sbatch -c 48 -p workers --wrap "hostname; cd /scratch && bash /home/guarracino/pggb-paper/scripts/gfa2evaluation.sh /lizardfs/guarracino/pggb-paper/graphs/scerevisiae8_p95.s1000.n8.k47.G4001-4507.SGDref/*.final.gfa SGDref /lizardfs/guarracino/pggb-paper/evaluation/scerevisiae8_p95.s1000.n8.k47.G4001-4507.SGDref $PATH_VCF_PREPROCESS $PATH_NUCMER_2_VCF 48"

sbatch -c 48 -p workers --wrap "hostname; cd /scratch && bash /home/guarracino/pggb-paper/scripts/gfa2evaluation.sh /lizardfs/guarracino/pggb-paper/graphs/scerevisiae8_p95.s5000.n8.k47.G4001-4507.SGDref/*.final.gfa SGDref /lizardfs/guarracino/pggb-paper/evaluation/scerevisiae8_p95.s5000.n8.k47.G4001-4507.SGDref $PATH_VCF_PREPROCESS $PATH_NUCMER_2_VCF 48"


sbatch -c 48 -p workers --wrap "hostname; cd /scratch && bash /home/guarracino/pggb-paper/scripts/gfa2evaluation.sh /lizardfs/guarracino/pggb-paper/graphs/ecoli100_p90.s5000.n100.k47.G4001-4507.BH100N_MG2017_3a53c37/*.final.gfa BH100N_MG2017_3a53c37 /lizardfs/guarracino/pggb-paper/evaluation/ecoli100_p90.s5000.n100.k47.G4001-4507.BH100N_MG2017_3a53c37/ $PATH_VCF_PREPROCESS $PATH_NUCMER_2_VCF 48"
```


Evaluation:

```shell
# Input
PATH_PGGB_VCF=/lizardfs/guarracino/pggb-paper/graphs/hsapiens90.chr6_p98.s5000.n90.k271.G4001-4507.chm13/hsapiens90.chr6.fa.gz.89d18c5.wfmash.paf.d2e0b63.1e53bae.smooth.final.grch38.vcf
LABEL=chr6_p98.s5000.n90.k271.G13117-13219

PATH_PGGB_VCF=/lizardfs/guarracino/pggb-paper/graphs/hsapiens90.chr6_p98.s5000.n90.k271.G13117-13219.chm13/hsapiens90.chr6.fa.gz.89d18c5.d2e0b63.59a83aa.smooth.final.grch38.vcf
LABEL=chr6_p98.s5000.n90.k271.G4001-4507

PATH_PGGB_VCF=/lizardfs/guarracino/pggb-paper/graphs/wgg88/wgg88.chr6.pan.fa.a2fb268.4030258.6a1ecc2.smooth.grch38.vcf
LABEL=wgg88.chr6

PATH_REF_FA=/lizardfs/guarracino/pggb-paper/references/hsapiens90.chr6.fa.gz.89d18c5.d2e0b63.59a83aa.smooth.final.grch38.fa
PATH_VCF_PREPROCESS=/home/guarracino/tools/pggb/scripts/vcf_preprocess.sh
DIR_REGIONS=/lizardfs/guarracino/HPRC/mini_dataset/union
DIR_TRUTH_VCF=/lizardfs/guarracino/HPRC/mini_dataset/truth
THREADS=48




mkdir -p $LABEL
cd $LABEL

echo "VCF file preprocessing for each sample"
grep '#CHROM' "$PATH_PGGB_VCF" -m 1 | cut -f 10- | tr '\t' '\n' | while read SAMPLE; do
    echo $SAMPLE

    bash "$PATH_VCF_PREPROCESS" \
        "$PATH_PGGB_VCF" \
        "$SAMPLE" \
        50 \
        "$PATH_REF_FA"
done

echo "Small variant evaluation"

PATH_REF_SDF="$PATH_REF_FA".sdf
if [[ ! -d $PATH_REF_SDF ]]; then
        rtg format -o "$PATH_REF_SDF" "$PATH_REF_FA"
fi

grep '#CHROM' "$PATH_PGGB_VCF" -m 1 | cut -f 10- | tr '\t' '\n' | while read SAMPLE; do
    echo $SAMPLE

    TRUTH_VCF=${DIR_TRUTH_VCF}/$SAMPLE.GRCh38_no_alt.deepvariant.prefix.vcf.gz

    if [[ -s $TRUTH_VCF ]]; then
        PATH_PGGB_SAMPLE_VCF="$PATH_PGGB_VCF"."$SAMPLE".max50.vcf.gz
        
        rtg vcfeval \
                -t $PATH_REF_SDF \
                -b $TRUTH_VCF \
                -c $PATH_PGGB_SAMPLE_VCF \
                -e <(zcat ${DIR_REGIONS}/GRCh38_notinalldifficultregions.bed.gz | sed 's/chr/grch38#1#chr/' ) \
                -T ${THREADS} \
                -o vcfeval/${SAMPLE}/easy

        rtg vcfeval \
                -t $PATH_REF_SDF \
                -b $TRUTH_VCF \
                -c $PATH_PGGB_SAMPLE_VCF \
                -e <(zcat ${DIR_REGIONS}/GRCh38_alldifficultregions.bed.gz | sed 's/chr/grch38#1#chr/') \
                -T ${THREADS} \
                -o vcfeval/${SAMPLE}/hard
    fi
done


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
```
