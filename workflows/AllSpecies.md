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
f=/lizardfs/guarracino/pggb-paper/sequences/hsapiens90.chr6.fa.gz
p=98
s=5000
n=90
k=271
G=13117,13219
ref=chm13
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c 48 -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -V $ref:# -t $t -o $out; mv $out $out_folder"

f=/lizardfs/guarracino/pggb-paper/sequences/hsapiens90.chr6.fa.gz
p=98
s=5000
n=90
k=271
G=4001,4507
ref=chm13
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c 48 -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -V $ref:# -t $t -o $out; mv $out $out_folder"

f=/lizardfs/guarracino/pggb-paper/sequences/hsapiens90.chr6.masked.fa.gz
p=98
s=5000
n=90
k=79
G=4001,4507
ref=chm13
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c 48 -w octopus09 --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O $O -G $G -V $ref:# -t $t -o $out; mv $out $out_folder"

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
s=1000
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




Evaluation:


```shell
vcfpreprocess=/home/guarracino/tools/pggb/scripts/vcf_preprocess.sh
nucmer2vcf=/home/guarracino/tools/pggb/scripts/nucmer2vcf.R


sbatch -c 48 -p workers --wrap "hostname; cd /scratch && bash /home/guarracino/pggb-paper/scripts/gfa2evaluation.sh /lizardfs/guarracino/pggb-paper/graphs/rat32.chr1_p98.s1000.n32.k271.G4001-4507.rn7/*.final.gfa rn7 /lizardfs/guarracino/pggb-paper/evaluation/rat32.chr1_p98.s1000.n32.k271.G4001-4507.rn7/ $vcfpreprocess $nucmer2vcf 48"

sbatch -c 48 -p workers --wrap "hostname; cd /scratch && bash /home/guarracino/pggb-paper/scripts/gfa2evaluation.sh /lizardfs/guarracino/pggb-paper/graphs/primates4.chr6_p95.s5000.n4.k47.G4001-4507.GRC38/*.final.gfa GRC38 /lizardfs/guarracino/pggb-paper/evaluation/primates4.chr6_p95.s5000.n4.k47.G4001-4507.GRC38/ $vcfpreprocess $nucmer2vcf 48"

sbatch -c 48 -p workers --wrap "hostname; cd /scratch && bash /home/guarracino/pggb-paper/scripts/gfa2evaluation.sh /lizardfs/guarracino/pggb-paper/graphs/primates4.chr7_p95.s5000.n4.k47.G4001-4507.GRC38/*.final.gfa GRC38 /lizardfs/guarracino/pggb-paper/evaluation/primates4.chr7_p95.s5000.n4.k47.G4001-4507.GRC38/ $vcfpreprocess $nucmer2vcf 48"

sbatch -c 48 -p workers --wrap "hostname; cd /scratch && bash /home/guarracino/pggb-paper/scripts/gfa2evaluation.sh /lizardfs/guarracino/pggb-paper/graphs/athaliana7.chr1_p95.s5000.n7.k47.G4001-4507.TAIR10/*.final.gfa TAIR10 /lizardfs/guarracino/pggb-paper/evaluation/athaliana7.chr1_p95.s5000.n7.k47.G4001-4507.TAIR10/ $vcfpreprocess $nucmer2vcf 48"

sbatch -c 48 -p workers --wrap "hostname; cd /scratch && bash /home/guarracino/pggb-paper/scripts/gfa2evaluation.sh /lizardfs/guarracino/pggb-paper/graphs/soy37.chr18_p95.s5000.n37.k47.G4001-4507.ZH13/*.final.gfa ZH13 /lizardfs/guarracino/pggb-paper/evaluation/soy37.chr18_p95.s5000.n37.k47.G4001-4507.ZH13 $vcfpreprocess $nucmer2vcf 48"

sbatch -c 48 -p workers --wrap "hostname; cd /scratch && bash /home/guarracino/pggb-paper/scripts/gfa2evaluation.sh /lizardfs/guarracino/pggb-paper/graphs/tomato23.chr2_p95.s5000.n23.k47.G4001-4507.SL5/*.final.gfa SL5 /lizardfs/guarracino/pggb-paper/evaluation/tomato23.chr2_p95.s5000.n23.k47.G4001-4507.SL5/ $vcfpreprocess $nucmer2vcf 48"


sbatch -c 48 -p workers --wrap "hostname; cd /scratch && bash /home/guarracino/pggb-paper/scripts/gfa2evaluation.sh /lizardfs/guarracino/pggb-paper/graphs/scerevisiae8_p95.s1000.n8.k47.G4001-4507.SGDref/*.final.gfa SGDref /lizardfs/guarracino/pggb-paper/evaluation/scerevisiae8_p95.s1000.n8.k47.G4001-4507.SGDref $vcfpreprocess $nucmer2vcf 48"

sbatch -c 48 -p workers --wrap "hostname; cd /scratch && bash /home/guarracino/pggb-paper/scripts/gfa2evaluation.sh /lizardfs/guarracino/pggb-paper/graphs/scerevisiae8_p95.s5000.n8.k47.G4001-4507.SGDref/*.final.gfa SGDref /lizardfs/guarracino/pggb-paper/evaluation/scerevisiae8_p95.s5000.n8.k47.G4001-4507.SGDref $vcfpreprocess $nucmer2vcf 48"


sbatch -c 48 -p workers --wrap "hostname; cd /scratch && bash /home/guarracino/pggb-paper/scripts/gfa2evaluation.sh /lizardfs/guarracino/pggb-paper/graphs/ecoli100_p90.s5000.n100.k47.G4001-4507.BH100N_MG2017_3a53c37/*.final.gfa BH100N_MG2017_3a53c37 /lizardfs/guarracino/pggb-paper/evaluation/ecoli100_p90.s5000.n100.k47.G4001-4507.BH100N_MG2017_3a53c37/ $vcfpreprocess $nucmer2vcf 48"
```
