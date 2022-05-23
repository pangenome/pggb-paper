## header ---------------------------------------------------------------------

library(data.table)

## settings -------------------------------------------------------------------

### fixed settings
argsVal <- commandArgs(trailingOnly = T)
dirBase <- argsVal[1]
dirOut <- file.path(dirBase, "os-vcf")
dir.create(path = dirOut, showWarnings = F, recursive = T)
inDir <- file.path(dirBase, "nuc-aln")

### Buff_noaln is [BUFF]: the distance from this SNP to the nearest mismatch
### (end of alignment, indel, SNP, etc) in the same alignment.
### Dist_seqend is [DIST]: the distance from this SNP
### to the nearest sequence end.
hdShowSnp <- c("Pos_a1", "Allele_a1", "Allele_a2", "Pos_a2",
               "Buff_noaln", "Dist_seqend", "Strand_a1", "Strand_a2",
               "Chrom_a1", "Chrom_a2")
hdShowSnpNoC <- c("Pos_a1", "Allele_a1", "Allele_a2", "Pos_a2",
                  "Buff_noaln", "Dist_seqend", "Nrep_a1", "Nrep_a2",
                  "Strand_a1", "Strand_a2", "Chrom_a1", "Chrom_a2")
allChr <- c("chrI", "chrII", "chrIII", "chrIV",
            "chrV", "chrVI", "chrVII", "chrVIII",
            "chrIX", "chrX", "chrXI", "chrXII",
            "chrXIII", "chrXIV", "chrXV", "chrXVI", "chrMT")

## clmnt ----------------------------------------------------------------------

stringInfo <- readLines(con = file.path(dirBase, "aux", "info-vcf.txt"),)

allFiles <- list.files(path = inDir, pattern = "var.txt$", full.names = T)

for (indF in allFiles) {
  ### set output file
  fileOut <- sub(pattern = "-var\\.txt", replacement = "-snps.vcf",
                 x = basename(indF))
  pathOut <- file.path(dirOut, fileOut)
  dtNucmer <- fread(indF, header = F, sep = "\t")
  colnames(dtNucmer) <- hdShowSnp
  ### filter out variants from non-unique alignments, and indels
  dtSNPs <- dtNucmer[Allele_a1 != "."
                     & Allele_a2 != "."]
  nR <- nrow(dtSNPs)
  dtVCF <- data.table(dtSNPs$Chrom_a1, dtSNPs$Pos_a1, rep(".", nR),
                      dtSNPs$Allele_a1, dtSNPs$Allele_a2, rep(60, nR),
                      rep(".", nR), # FILTER
                      rep(".", nR), # INFO
                      rep("GT:FRMR:FRMQ:CHRQ:STARTQ:ENDQ", nR), # FORMAT
                      paste("1", dtSNPs$Strand_a1, dtSNPs$Strand_a2,
                            dtSNPs$Chrom_a2, dtSNPs$Pos_a2,
                            dtSNPs$Pos_a2, sep = ":"))
  ### sorting by chromosome and position
  dtVCF <- dtVCF[order(match(dtVCF$V1, allChr), dtVCF$V2), ]

  ## make the header ----------------------------------------------------------

  headerVcf <- c("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
  nameSample <- basename(indF)
  idSample <- paste(unlist(strsplit(nameSample, split = "-"))[4:5],
                    collapse = "-")
  string1 <- "##fileformat=VCFv4.3"
  string2 <- paste0("##fileDate=", date())
  string3 <- paste0("##source=nucmer v4.0.0")
  stringH <- paste0(headerVcf, idSample)
  nameRef <- paste(unlist(strsplit(nameSample, split = "-"))[1:2],
                   collapse = "-")
  pathRef <- list.files(path = file.path(dirBase, "genomes", "reference"),
                        pattern = paste0(nameRef, ".*fa$"),
                        full.names = T, recursive = F)
  stringRef <- paste0("##reference=", pathRef)
  pathRefFai <- paste0(pathRef, ".fai")

  dtFai <- fread(pathRefFai, header = F,
                 sep = "\t")[, c(1, 2)]
  colnames(dtFai) <- c("Chrom_id", "Len_bp")
  stringContig <- paste0("##contig=<ID=",
                         dtFai$Chrom_id, ",length=", dtFai$Len_bp, ">")

  cat(string1, string2, string3,
      stringRef,
      stringContig,
      stringInfo,
      stringH,
      file = pathOut, sep = "\n")
  ### append the vcf to the corresponding header (with the correct sample name)
  fwrite(file = pathOut, append = T, x = dtVCF, quote = F,
              sep = "\t", row.names = F, col.names = F)
}
