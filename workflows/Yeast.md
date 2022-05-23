# Yeast

## Preamble

In this tutorial we will use the SGD reference genome and 7 yeast genomes from 7 strains of the *S. cerevisiae* to build a variation graph and call variants from it. The aforementioned 7 strains are described in Yue *et al.* [Yue J, ..., & Liti G. Contrasting evolutionary genome dynamics between domesticated and wild yeasts. NAT GENET. 2017.]. They represent different populations of homozygous diploids which were sequenced with PacBio and Illumina platforms and *de novo* assembled. The yeast reference genome is described here: https://www.yeastgenome.org.

The biallelic SNVs called from the variation graph will be compared to those obtained from the pairwise alignments of each the 7 genomes against the SGD reference. The latter will be used as the ground truth to calculate precision, recall, as well as the F1-score. Moreover, the ground truth will be stratified in two types of variants: those lying is "easy" regions and those lying in "hard regions". The latter include genomic regions which are repeated or features, such as Tys (yeast transposable elements), which are known to be repeated.

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

Here, we align the SGD genome against any other assembly. Before running the following bash lines, "cd" to the main folder of the yeast repository. This chunk produces the *-coords.txt and *-var.txt for each of the genomes aligned against SGD. The former is the file  with the detected variants while the latter is the file reporting the alignments. Both files can be found in the "nuc-aln" folder which has been created in the main folder of the repository.

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

Now, we  have to convert the output of nucmer to a vcf. For this purpose we use a custom R script (nucmer-vcf.R) which is in the "scripts" folder. To run it you can use the following command line from the main folder.

```
Rscript scripts/nucmer-vcf.R $(pwd)
```

The single-sample vcf files can be converted in a multi-sample vcf file.

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

### prepare header for the next chunk
bgzip -b -c "${dir_out}/multis-snps.vcf.gz" | grep "^#" \
> "${dir_out}/multis-snps-genfix.vcf"

### back to the main folder
cd "${dir_multis}"
```

Now, we have to check whether the genotypes which are missing (e.g. any entry reporting "." in the fields of columns 10-16 of the multi-sample vcf file) are due to the presence of a reference allele or to a missing alignment. In the first case we replace them with a reference genotype (namely "0") otherwise we keep the missing genotype (".") record. For this purpose, we use an R script stored in the "script" folder.

```
Rscript scripts/check-non-alt.R $(pwd)
```

## Graph construction and variant calls
