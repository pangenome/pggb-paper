# Yeast

## Preamble

In this tutorial we will use the SGD reference genome and 7 yeast genomes from 7 strains of the *S. cerevisiae* to build a variation graph and call variants from it. The aforementioned 7 strains are described in Yue *et al.* [Yue J, ..., & Liti G. Contrasting evolutionary genome dynamics between domesticated and wild yeasts. NAT GENET. 2017.]. They represent different populations of homozygous diploids which were sequenced with PacBio and Illumina platforms and *de novo* assembled. The yeast reference genome is described here: https://www.yeastgenome.org.

The variants called from the variation graph will be compared to those obtained from the pairwise alignments of each the 7 genomes against the SGD reference. The latter will be used as the ground truth to calculate precision, recall, as well as the F1-score. Moreover, the ground truth will be stratified in two types of variants: those lying is "easy" regions and those lying in "hard regions". The latter include genomic regions which are repeated or features, such as Tys (yeast transposable elements), which are known to be repeated.

All the input sequences are provided in the "genomes" folder while the corresponding annotations are stored in the "annotations" folder. All the scripts needed to run the analysis are stored in the "scripts" folder.

## Data

The structure of the "genomes" folder is:
```
genomes/
├── collinear
├── rearranged
└── reference
```
The genome are organised using *a priori* knowledge (as reported in Yue *et al.*). The files are already well-formatted for the following scripts. Yeast's chromosomes are encoded with Roman numerals, so e.g. chromosome ten is encoded as "chrX". Nobody likes it but do not ask a yeast geneticist to change it.

## Dependencies

The following command-line tools and R packages are required:

- bcftools
- bgzip
- tabix
- mummer
- wait -n (may be missing in macOS)
- bedtools
- RepeatMasker
- rtg
- data.table (R package)

## Graph construction and variant calls
