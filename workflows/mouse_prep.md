Prepare references:

```shell
# Take chromosomes from 1-22, X, Y, and MT

mkdir -p /lizardfs/guarracino/pggb-paper/references

cd /lizardfs/erikg/mouse/assemblies
cat \
    <( samtools faidx Mus_musculus.GRCm38.68.dna.toplevel.fa.gz $(head -n 22 Mus_musculus.GRCm38.68.dna.toplevel.fa.gz.fai | cut -f 1) ) \
    <( samtools faidx Mus_musculus.GRCm39.dna.toplevel.fa.gz $(head -n 22 Mus_musculus.GRCm39.dna.toplevel.fa.gz.fai | cut -f 1) ) | \
    bgzip -c -@ 48 > /lizardfs/guarracino/pggb-paper/references/GRCm38+GRCm39.fa.gz
samtools faidx /lizardfs/guarracino/pggb-paper/references/GRCm38+GRCm39.fa.gz
```

Map contigs against the references:

```shell
mkdir -p /lizardfs/guarracino/pggb-paper/assemblies/mouse_partitioning/
cd /lizardfs/guarracino/pggb-paper/assemblies/mouse_partitioning/

REFERENCES_FASTA=/lizardfs/guarracino/pggb-paper/references/GRCm38+GRCm39.fa.gz

(ls /lizardfs/erikg/mouse/assemblies/*fasta.gz) | while read FASTA; do
  HAPLOTYPE=$(basename $FASTA .fa.gz | cut -f 1,2 -d '_');
  echo $HAPLOTYPE
  
  PAF=/lizardfs/guarracino/pggb-paper/assemblies/mouse_partitioning/$HAPLOTYPE.vs.refs.paf
  sbatch -p workers -c 12 --wrap "wfmash -t 12 -m -s 5k -p 90 -H 0.001 $REFERENCES_FASTA $FASTA > $PAF"
done
```

Get the best target for each contig:


```shell
(ls /lizardfs/erikg/mouse/assemblies/*fasta.gz) | while read FASTA; do
  HAPLOTYPE=$(basename $FASTA .fa.gz | cut -f 1,2 -d '_');
  echo $HAPLOTYPE
  
  PAF=/lizardfs/guarracino/pggb-paper/assemblies/mouse_partitioning/$HAPLOTYPE.vs.refs.paf

  # We remove the reference prefixes because we are interested in the chromosome
  python3 /lizardfs/guarracino/chromosome_communities/scripts/partition_by_weighted_mappings.py <(sed 's/GRCm38#1#//g;s/GRCm39#1#//g' $PAF) > $HAPLOTYPE.vs.refs.partitions.tsv
done
```

Take chromosome 19:

```shell
cd /lizardfs/guarracino/pggb-paper/sequences/

REFERENCES_FASTA=/lizardfs/guarracino/pggb-paper/references/GRCm38+GRCm39.fa.gz

samtools faidx $REFERENCES_FASTA $(grep chr19 $REFERENCES_FASTA.fai | cut -f 1) > mouse17.chr19.fa

(ls /lizardfs/erikg/mouse/assemblies/*fasta.gz) | while read FASTA; do
  HAPLOTYPE=$(basename $FASTA .fa.gz | cut -f 1,2 -d '_');
  echo $HAPLOTYPE
  
  PARTITIONS=/lizardfs/guarracino/pggb-paper/assemblies/mouse_partitioning/$HAPLOTYPE.vs.refs.partitions.tsv
  
  samtools faidx $FASTA $( awk '$2 == "chr19"' < $PARTITIONS | cut -f 1 ) >> mouse17.chr19.fa
done

bgzip -@ 48 mouse17.chr19.fa
samtools faidx mouse17.chr19.fa.gz
```
