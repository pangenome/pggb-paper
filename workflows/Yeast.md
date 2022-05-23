# Yeast

## Preamble

In this tutorial we will use 7 yeast genomes from 7 strains of the *S. cerevisiae* to build a variation graph and call variants from it. The strains are described in Yue *et al.* [Yue J, ..., & Liti G. Contrasting evolutionary genome dynamics between domesticated and wild yeasts. NAT GENET. 2017.]. They represent different populations of homozygous diploids which were sequenced with PacBio and Illumina platforms.

The variants called from the variation graph will be compared to those called comparing the genomes by pairwise alignment (which will be used as the ground truth to calculate precision, recall, as well as the F1-score). Moreover, the ground truth will be stratified in two types of variants: those lying is "easy" regions and those lying in "hard regions". The latter include genomic regions which are repeated or features, such as Tys (yeast transposable elements), which are known to be repeated.

All the input data are provided in the "data" folder. All the scripts need to run the analysis are stored in the "scripts" folder.

## Data

The tree of the "data" folder is:
```
data
├── annotations
│   ├── mitochondrial
│   └── nuclear
└── genomes
    ├── haploid-collapsed
    └── mitochondrial
```

The annotation are provided as GFF files. Mitochondrial and nuclear data are separated for convenience. The same holds for genomes. The latter are provided as *.fa.gz files. Indexes (obtained from "samtools faidx" of the fasta files) are also included. The files are already well-formatted for the following scripts. Yeast's chromosomes are encoded with Roman numerals, so e.g. chromosome ten is encoded as "chrX". Nobody likes it but do not ask a yeast geneticist to change it.

## Graph construction and variant calls
