## header ---------------------------------------------------------------------

library(data.table)

## settings -------------------------------------------------------------------

### fixed settings
argsVal <- commandArgs(trailingOnly = T)
dirBase <- argsVal[1]
dirOut <- file.path(dirBase, "ms-vcf")
dir.create(path = dirOut, showWarnings = F, recursive = T)
dirIn <- file.path(dirBase, "ms-vcf")
dirAln <- file.path(dirBase, "nuc-aln")
fileOut <- "multis-snps-genfix.vcf"
hdAlnCoords <- c("Start_a1", "End_a1", "Start_a2", "End_a2", "Len_a1", "Len_a2",
                 "Perc_iden", "Chrom_a1", "Chrom_a2")

## clmnt ----------------------------------------------------------------------

fileIn <- list.files(path = dirIn, pattern = "vcf.gz$", full.names = T)
dtInVcf <- fread(fileIn)
colnames(dtInVcf)[1] <- "CHROM"

### remove non-biallelic loci
indNonBiallelic <- grep(pattern = ",", x = dtInVcf$ALT)
nNonBiallelic <- length(indNonBiallelic)
cat("Found ", nNonBiallelic, " non-biallelic SNPs out of ",
    nrow(dtInVcf), "\n", sep = "")
dtBiallelicVcf <- dtInVcf[!indNonBiallelic]

nColBiallelic <- ncol(dtBiallelicVcf)

sampNames <- colnames(dtBiallelicVcf)[10:nColBiallelic]

### count FORMAT tags and set the new string
nTags <- length(unlist(strsplit(dtBiallelicVcf$FORMAT[1], split = ":")))
newTag <- paste(c("0", rep(".", c(nTags - 1))), collapse = ":")

for (indS in sampNames) {
  ### indexes of column indS to be checked for the REF allele
  indC <- grep(pattern = "^\\.", x = dtBiallelicVcf[[indS]])
  cat("There are ", length(indC), " alleles to be checked for sample ", indS,
      "\n", sep = "")
  ### reading the alignment
  strAln <- paste0("SGDref-hc", "-vs-", indS, "-coords.txt")
  fileAln <- list.files(dirAln, pattern = strAln, full.names = T)
  if (!file.exists(fileAln)) {
    cat("File ", fileAln, "does not exist. Check the IDs of the genomes!",
        "\n", sep = "")
    stop()
  }
  dtAln <- fread(fileAln, sep = "\t")
  colnames(dtAln) <- hdAlnCoords

  ### subset vcf and format loci for the comparison
  dtLoci <- data.table(Chr_ref = dtBiallelicVcf$CHROM[indC],
                       Start_ref = dtBiallelicVcf$POS[indC],
                       End_ref = dtBiallelicVcf$POS[indC],
                       Allele_ref = dtBiallelicVcf$REF[indC])
  setkey(dtLoci, Chr_ref, Start_ref, End_ref)
  ### format alignment for the comparison
  setkey(dtAln, Chrom_a1, Start_a1, End_a1)

  ### make the comparison using keys from the vcf file
  dtOverlap <- foverlaps(dtLoci, dtAln, nomatch = NA, mult = "all",
                          by.x = c("Chr_ref", "Start_ref", "End_ref"))
  ### reduce the overlap table removing loci that overlap multiple alignments
  boDup <- duplicated(dtOverlap, by = c("Chr_ref", "Start_ref", "End_ref"))
  dtOverlapReduced <- dtOverlap[!boDup, ]
  ### indexes of SNPs with multiple alignments in the reduced table
  dtMultiAln <- dtOverlap[boDup, ]
  strOverlapReduced <- paste(dtOverlapReduced$Chr_ref,
                             dtOverlapReduced$Start_ref,
                             sep = "_")
  strMultiAln <- paste(dtMultiAln$Chr_ref,
                       dtMultiAln$Start_ref,
                       sep = "_")
  indMultiAln <- match(strMultiAln, strOverlapReduced)
  ### boolean of SNPs within one or more alignments (TRUE)
  boAln <- !is.na(dtOverlapReduced$Perc_iden)
  ### boolean of SNPs within only one alignment (TRUE)
  boAln[indMultiAln] <- F

  ### update the GT tag
  nCol <- which(colnames(dtBiallelicVcf) == indS)
  dtBiallelicVcf[indC[boAln], nCol] <- newTag
}

fwrite(x = dtBiallelicVcf, file = file.path(dirOut, fileOut), append = T,
       quote = F, sep = "\t", row.names = F, col.names = F)
