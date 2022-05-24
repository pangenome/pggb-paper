# Yeast

## Preamble

In this tutorial we will use the *S. cerevisiae* SGD reference genome and 7 yeast assemblies from 7 strains of the *S. cerevisiae* to build a variation graph and call variants from it. The aforementioned 7 strains are described in Yue *et al.* [Yue J, ..., & Liti G. Contrasting evolutionary genome dynamics between domesticated and wild yeasts. NAT GENET. 2017.]. They represent different populations of homozygous diploids which were sequenced with PacBio and Illumina platforms and *de novo* assembled. The yeast reference genome is described here: https://www.yeastgenome.org.

The bi-allelic SNVs called from the variation graph will be compared to those obtained from the pairwise alignments of each the 7 genomes against the SGD reference. The latter will be used as the ground truth to calculate precision, recall, as well as the F1-score. Moreover, the ground truth will be stratified in two types of variants: those lying is "easy" regions and those lying in "hard regions". The latter include genomic regions which are repeated or features, such as Tys (yeast transposable elements), which are known to be repeated.

All the input sequences are provided in the "genomes" folder while the corresponding annotations are stored in the "annotations" folder. All the scripts needed to run the analysis are stored in the "scripts" folder. Make sure to "cd" to the main folder of the repository before running the command lines reported below.

## Data

The structure of the "genomes" folder is:

```
genomes
├── collinear
├── rearranged
└── reference
```

The genomes are organised using *a priori* knowledge (as reported in Yue *et al.*). The files are already well-formatted for the following scripts. Yeast's chromosomes are encoded with Roman numerals, so e.g. chromosome ten is encoded as "chrX". Nobody likes it but do not ask a yeast geneticist to change it.

## Dependencies

The following command-line tools and R packages are required:

- pggb
- bcftools
- bgzip
- tabix
- mummer
- wait -n (may be missing in macOS)
- bedtools
- RepeatMasker
- rtg
- data.table (R package)

## Clone the yeast repository

First, we have to clone the yeast repository.

```
git clone -b yeast https://github.com/pangenome/pggb-paper.git
```

## Variant call from the linear genomes

Here, we align the SGD genome against any other assembly. Before running the following bash lines, "cd" to the main folder of the yeast repository. This chunk of code produces the *-coords.txt and *-var.txt for each of the genomes aligned against SGD. The former is the file  with the detected variants while the latter is the file reporting the alignments. Both files can be found in the "nuc-aln" folder which has been created in the main folder of the repository.

```
## settings -------------------------------------------------------------------

### base folder
dir_aln=$(pwd)
### parallel runs
pll_runs=4
### suppress long messages from background processes
set +m

## clmnt ----------------------------------------------------------------------

### clean and make output folder
out_dir="${dir_aln}/nuc-aln"
if [[ -d "${out_dir}" ]]; then
  rm -rf "${out_dir}"
fi
mkdir -p "${out_dir}"

### reference file
ref_dir="${dir_aln}/genomes/reference"
ref_seq=$(find "${ref_dir}" -name "*fa")
ref_name=$(basename "${ref_seq}" | sed 's|-genome\.fa||')

## collinear assemblies -------------------------------------------------------

in_dir="${dir_aln}/genomes/collinear"
cd "${in_dir}"
if [[ -d "tmp" ]]; then
  rm -rf "tmp"
fi
mkdir "tmp"

### alignment against SGD
seq_arr=( $(find . -name "*fa") )
seq_dim=$(echo "${#seq_arr[@]}")
ref_chroms=$(grep ">chr" "${ref_seq}" | sed 's|>||g')

for (( ind_i=0; ind_i<seq_dim; ind_i++ )); do
  name_i=$(echo "${seq_arr[ind_i]}" | sed 's|-genome\.fa||' | sed 's|^..||')
  out_prefix="${ref_name}-vs-${name_i}"

  ### nucmer chromosome-wise
  for ind_c in ${ref_chroms}; do
    ((cnt_p++))
    (
    ### check if the query has the contig
    check_contig=$(grep -w "${ind_c}" "${seq_arr[ind_i]}.fai")
    if [[ -z "${check_contig}" ]]; then
      echo "[alignment with nucmer]  ${ind_c} not found in ${seq_arr[ind_i]}.fai"
      continue
    fi

    ### extract chromosome sequence
    samtools faidx ${ref_seq} "${ind_c}" > "tmp/${ref_name}-${ind_c}.fa"
    samtools faidx ${seq_arr[ind_i]} "${ind_c}" > "tmp/${name_i}-${ind_c}.fa"

    ### go nucmer, go!
    nucmer --prefix="${out_dir}/${out_prefix}-${ind_c}" \
    "tmp/${ref_name}-${ind_c}.fa" "tmp/${name_i}-${ind_c}.fa"

    show-snps -THC "${out_dir}/${out_prefix}-${ind_c}.delta" \
    > "${out_dir}/${out_prefix}-${ind_c}-var.txt"

    show-coords -TH "${out_dir}/${out_prefix}-${ind_c}.delta" \
    > "${out_dir}/${out_prefix}-${ind_c}-coords.txt"
    ) &
  if (( cnt_p % pll_runs == 0 )); then
    wait -n
    cnt_p=$(( pll_runs - 1 ))
  fi
  done
done

### wait for all chromosomes to be processed
wait

### clean tmp folder
rm -rf "tmp"

### merge the single-chromosome data in a whole-genome file
for ind_f in $(find . -name "*fai"); do
  cont_sort=$(cut -f 1 "${ind_f}")
  strain_id=$(echo "${ind_f}" | cut -d "-" -f 1,2 | sed 's|^..||')
  out_var="${ref_name}-vs-${strain_id}-var.txt"
  out_coords="${ref_name}-vs-${strain_id}-coords.txt"
  ### cleaning
  if [[ -f "${out_dir}/${out_var}" ]]; then
    rm -f "${out_dir}/${out_var}"
  fi
  if [[ -f "${out_dir}/${out_coords}" ]]; then
    rm -f "${out_dir}/${out_coords}"
  fi
  ### merging
  for ind_c in ${cont_sort}; do
    cat "${out_dir}/${ref_name}-vs-${strain_id}-${ind_c}-var.txt" \
    >> "${out_dir}/${out_var}"
    cat "${out_dir}/${ref_name}-vs-${strain_id}-${ind_c}-coords.txt" \
    >> "${out_dir}/${out_coords}"
    ### clean single-chromosome data
    rm -f "${out_dir}/${ref_name}-vs-${strain_id}-${ind_c}-var.txt"
    rm -f "${out_dir}/${ref_name}-vs-${strain_id}-${ind_c}-coords.txt"
  done
done

## rearranged assemblies ------------------------------------------------------

in_dir="${dir_aln}/genomes/rearranged"
cd "${in_dir}"

### alignment against SGD
seq_arr=( $(find . -name "*fa") )
seq_dim=$(echo "${#seq_arr[@]}")

for (( ind_i=0; ind_i<seq_dim; ind_i++ )); do
  ((cnt_p++))
  (
  name_i=$(echo "${seq_arr[ind_i]}" | sed 's|-genome\.fa||' | sed 's|^..||')
  out_prefix="${ref_name}-vs-${name_i}"

  ### go nucmer, go!
  nucmer --prefix="${out_dir}/${out_prefix}-wg" \
  "${ref_seq}" ${seq_arr[ind_i]}

  # delta-filter -u 100 "${out_dir}/${out_prefix}-wg.delta" \
  # > "${out_dir}/${out_prefix}-wg-flt.delta"

  show-snps -THC "${out_dir}/${out_prefix}-wg.delta" \
  > "${out_dir}/${out_prefix}-var.txt"

  show-coords -TH "${out_dir}/${out_prefix}-wg.delta" \
  > "${out_dir}/${out_prefix}-coords.txt"

  ### delta file
  # rm -f "${out_dir}/${out_prefix}.delta"
  ) &
  if (( cnt_p % pll_runs == 0 )); then
    wait -n
    cnt_p=$(( pll_runs - 1 ))
  fi
done

### clean all delta files (both for collinear and rearranged genomes)
find "${out_dir}" -name "*delta" | xargs rm

### back to the main folder
cd "${dir_aln}"
```

Now, we  have to convert the output of nucmer to a vcf. For this purpose we use a custom R script (nucmer-vcf.R) which is in the "scripts" folder. To run it you can use the following command line from the main folder. The results are stored in the "os--vcf" folder.

```
Rscript scripts/nucmer-vcf.R $(pwd)
```

The single-sample vcf files can be converted in a multi-sample vcf file. It is stored in the "ms-folder".

```
## settings -------------------------------------------------------------------

### threads
n_threads=4

### base dir
dir_multis=$(pwd)

### input and output folders
dir_input="${dir_multis}/os-vcf"
dir_out="${dir_multis}/ms-vcf"
if [[ -d "${dir_out}" ]]; then
  rm -rf "${dir_out}"
fi
mkdir -p "${dir_out}"

## clmnt ----------------------------------------------------------------------

cd "${dir_input}"
rm -f *tbi

### compress with bgzip all the single-sample vcf files and index with tabix
for vcf_file in $(find . -name "*vcf"); do
  bgzip -c -f "${vcf_file}" > "${vcf_file}.gz"
  tabix "${vcf_file}.gz"
done

### merge to one multisample vcf with bcftools
all_files=( $(find . -name "*vcf.gz") )
bcftools merge --output-type z --threads "${n_threads}" ${all_files[@]} \
> "${dir_out}/multis-snps.vcf.gz"

### prepare header for the next chunk of code
bgzip -b -c "${dir_out}/multis-snps.vcf.gz" | grep "^#" \
> "${dir_out}/multis-snps-genfix.vcf"

### back to the main folder
cd "${dir_multis}"
```

Now, we have to check whether the genotypes which are missing (e.g. any entry reporting "." in the fields of columns 10-16 of the multi-sample vcf file) are due to the presence of a reference allele or to a missing alignment. In the first case we replace them with a reference genotype (namely "0") otherwise we keep the missing genotype (".") record. For this purpose, we use an R script stored in the "script" folder. The output vcf file (multis-snps-genfix.vcf) is stored in the "ms-vcf" folder.

```
Rscript scripts/check-non-alt.R $(pwd)
```
## Variants stratification

Calling variants in some specific regions of a genome can be particularly difficult due to the presence of repeated sequences. Here, we use RepeatMasker to make a bed files with easy and hard regions. Plus, we add any region corresponding to the following annotations in the gff file of the SGD reference genome: "CEN" (centromeres), "X_element" and "Y_prime_element" (X and Y' elements respectively), "TY" (yeast transposable elements), and "soloLTR" (solo long terminal repeats).

```
## settings -------------------------------------------------------------------

### base dir
dir_base=$(pwd)

### stratification folders
dir_out="${dir_base}/rep-mask"

### other variables
ref_id="SGDref-hc"
ref_gz="${ref_id}-genome.fa.gz"
ref_fa="${ref_id}-genome.fa"
ref_path="${dir_base}/genomes/reference/${ref_fa}"
feat_gff="${dir_base}/annotations/${ref_id}-features.gff"

### clean working folder
if [[ -d "${dir_out}" ]]; then
  rm -rf "${dir_out}"
fi
mkdir -p "${dir_out}"

## clmnt ----------------------------------------------------------------------

# to run RepeatMasker (default ReapeatMasker), with also the MT genome
RepeatMasker -species fungi -xsmall -gff "${ref_path}" -dir "${dir_out}"

### format the gff (version 2)
cd "${dir_out}"
grep -v "^#" "SGDref-hc-genome.fa.out.gff" | \
sed 's|RepeatMasker|SGDref|g' | \
sed 's|Target |Target=|g' | sed 's|"||g' \
> "reg-mask.gff"

### how many bp were masked by RepeatMasker
echo "Number of bp masked by RepeatMasker:"
grep -v "^#" "reg-mask.gff" \
| cut -f 4,5 \
| awk 'BEGIN {FS="\t"} { sum += $2-$1 } END {print sum}'

### get hard features from the gff file of SGD (both nuclear and MT)
### things to be extracted from the annotation of SGD
### and marked as hard regions:
### CEN
### X_element
### Y_prime_element
### TY
### [S/s]olo_LTR
grep "^#" "${feat_gff}" > "feat-hard.gff"
str_feat="CEN X_element Y_prime_element TY oloLTR"
for ind_s in ${str_feat}; do
  grep "${ind_s}" ${feat_gff} >> "feat-hard.gff"
done

### how many bp were masked by feature-selection
echo "Number of bp found in hard-feature regions:"
grep -v "^#" "feat-hard.gff" \
| cut -f 4,5 \
| awk 'BEGIN {FS="\t"} { sum += $2-$1 } END {print sum}'

### compare the two sets: overlap in bp
echo "Overlap between masked and hard-feature regions:"
bedtools intersect -a "feat-hard.gff" \
-b "reg-mask.gff" | \
awk 'BEGIN {FS="\t"} { sum += $5-$4 } END { print sum }'

### merge the two gff
cat "feat-hard.gff" "reg-mask.gff" \
> "reg-hard.gff"

### sort
all_chr="chrI chrII chrIII chrIV chrV chrVI chrVII chrVIII \
chrIX chrX chrXI chrXII chrXIII chrXIV chrXV chrXVI chrMT"
grep "^#" "reg-hard.gff" > "reg-hard-srt.gff"
for ind_c in ${all_chr}; do
  grep -w "${ind_c}" "reg-hard.gff" | \
  sort -k 4,4n >> "reg-hard-srt.gff"
done

### convert the gff to bed and merge
### (i.e. remove overlapping features) the intervals
cut -f 1,4,5 "reg-hard-srt.gff" | bedtools merge > "reg-hard-srt-mrg.bed"

### count the bp of the final bed (no overlapping feature is left)
echo "Number of bp in hard regions \
(no overlap between masked regions and hard-feature regions):"
awk 'BEGIN {FS="\t"} {sum += $3-$2} END {print sum}' "reg-hard-srt-mrg.bed"

### now we calculate the easy regions
all_reg_bed="all-reg.bed"
awk 'BEGIN {FS="\t"} {OFS="\t"} {print $1,1,$2}' "${ref_path}.fai" \
> "${all_reg_bed}"
bedtools subtract -a "${all_reg_bed}" -b "reg-hard-srt-mrg.bed" \
> "reg-easy-srt-mrg.bed"
```

## Graph construction and variant calls

## Calculation of precision, recall (sensitivity), and  the F-score

Only "callable" positions are taken into account.
