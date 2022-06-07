### All species

The (partial) datasets are in `/lizardfs/guarracino/pggb-paper/sequences`.

`pggb` command lines:

```shell
out_folder=/lizardfs/guarracino/pggb-paper/graphs

# Fixed parameters
t=40
POA=asm5
O=0.03

# Homo sapiens
f=/lizardfs/guarracino/pggb-paper/sequences/hsapiens90.chr6.fa.gz
p=98
s=5000
n=90
k=311
G=13117,13219
ref=chm13
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c 48 -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O 0$O -G $G -V $ref:# -t $t -o $out; mv $out $out_folder"

# Mus musculus
f=/lizardfs/guarracino/pggb-paper/sequences/mouse149.chr19.fa.gz
p=98
s=5000
n=149
k=229
G=4001,4507
ref=REF
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c 48 -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O 0$O -G $G -V $ref:# -t $t -o $out; mv $out $out_folder"

# Rat
f=/lizardfs/guarracino/pggb-paper/sequences/rat32.chr1.fa.gz
p=98
s=5000
n=32
k=229
G=4001,4507
ref=rn7
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c 48 -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O 0$O -G $G -V $ref:# -t $t -o $out; mv $out $out_folder"

# Primates
f=/lizardfs/guarracino/pggb-paper/primates4.chr7.fa.gz
p=95
s=5000
n=4
k=47
G=4001,4507
ref=GRC38
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c 48 -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O 0$O -G $G -V $ref:# -t $t -o $out; mv $out $out_folder"

# Arabidopsis thaliana
f=/lizardfs/guarracino/pggb-paper/athaliana7.chr1.fa.gz
p=95
s=5000
n=7
k=47
G=4001,4507
ref=TAIR10
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c 48 -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O 0$O -G $G -V $ref:# -t $t -o $out; mv $out $out_folder"

# Soy
f=/lizardfs/guarracino/pggb-paper/soy37.chr18.fa.gz
p=95
s=5000
n=37
k=47
G=4001,4507
ref=ZH13
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c 48 -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O 0$O -G $G -V $ref:# -t $t -o $out; mv $out $out_folder"

# Tomato
f=/lizardfs/guarracino/pggb-paper/tomato23.chr2.fa.gz
p=95
s=5000
n=23
k=47
G=4001,4507
ref=SL5
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c 48 -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O 0$O -G $G -V $ref:# -t $t -o $out; mv $out $out_folder"

# Yeast
f=/lizardfs/guarracino/pggb-paper/sequences/scerevisiae8.fa.gz
p=95
s=5000
n=8
k=47
G=4001,4507
ref=SGDref
out=$(basename "$f" .fa.gz)_p$p.s$s.n$n.k$k.G$(echo $G | tr ',' '-').$ref
sbatch -c 48 -p 386mem --wrap "hostname; cd /scratch; pggb -i $f -p $p -s $s -n $n -k $k -P $POA -O 0$O -G $G -V $ref:# -t $t -o $out; mv $out $out_folder"
```


Evaluation:


```shell
bash /home/guarracino/pggb-paper/scripts/gfa2evaluation.sh \
    $out_folder/.../...smooth.final.gfa \
    chm13 \
    $ref \
    /home/guarracino/tools/pggb/scripts/vcf_preprocess.sh \
    /home/guarracino/tools/pggb/scripts/nucmer2vcf.R \
    48
```
