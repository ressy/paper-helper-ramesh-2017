source("scripts/analysis.R")

# Go ----------------------------------------------------------------------


genbank <- read.csv("converted/all.csv", stringsAsFactors = FALSE)
genbank_alleles <- parse_genbank_genes(genbank)
sonar_alleles <- load_sonar_alleles()

# SONAR generally has fewer entries than what the paper described, but not for
# all (for example IGKV), and it includes ones labeled "partial" that the paper
# seems to exclude from the final counts.
# segments <- c(paste0("IGH", c("V", "D", "J")), "IGLV", "IGLJ", "IGKV", "IGKJ")
# table(subset(genbank_alleles, Segment %in% segments)$Segment)[segments]
# table(subset(genbank_alleles, Segment %in% segments & ! AccessionDescriptionPartial)$Segment)[segments]
# table(subset(sonar_alleles, Segment %in% segments)$Segment)[segments]

paper <- load_paper()
metadata <- list(
  genes = parse_genes(paper),
  scaffolds = parse_scaffolds(paper)
)

genbank_alleles <- merge_metadata(genbank_alleles, metadata)


# Checks ------------------------------------------------------------------


.genes_main <- function(x) {
  unique(subset(
    x,
    Dataset == "WGS" &
      GeneCategory == "Functional" &
      ScaffoldLocusGroup %in% c("main", "sister"))$Gene)
}

.alleles_main <- function(x) {
  subset(
    x,
    GeneCategory == "Functional" & 
      ! ScaffoldLocusGroup %in% "unknown")$Allele
}

.genes_unknown <- function(x) {
  unique(subset(
    x,
    Dataset == "WGS" &
      GeneCategory == "Functional" &
      ScaffoldLocusGroup %in% c(NA, "unknown"))$Gene)
}

.alleles_unknown <- function(x, genes_unknown) {
  subset(x, Gene %in% genes_unknown)$Allele
}

.genes_other <- function(x, genes_main, genes_unknown) {
  genes_other <- unique(subset(
    x,
    Dataset == "Targeted" &
      GeneCategory == "Functional")$Gene)
  genes_other <- genes_other[! genes_other %in% genes_main]
  genes_other <- genes_other[! genes_other %in% genes_unknown]
  genes_other
}

# Figure 4 shows five functional IGHV2 genes, but table 1 shows 6 (one with
# unknown scaffold).  Why?  I see no functional genes with unknown scaffold
# placement.
# Figure 4 shows 39 functional IGHV3 genes, but table 1 shows only 38.  I see
# all 39.
# TODO make the allele counts make sense
make_table1 <- function(genbank_alleles) {
  chunk <- subset(genbank_alleles, Segment == "IGHV")
  result <- do.call(rbind, lapply(split(chunk, chunk$Family), function(chunk_fam) {
    genes_main <- .genes_main(chunk_fam)
    genes_unknown <- .genes_unknown(chunk_fam)
    alleles_unknown <- .alleles_unknown(chunk_fam, genes_unknown)
    genes_other <- .genes_other(chunk_fam, genes_main, genes_unknown)
    alleles_main <- .alleles_main(chunk_fam)
    val1 <- c(length(genes_main), length(genes_other))
    val1 <- val1[val1 != 0]
    val1 <- paste0(paste(val1, collapse=" + "), " (", length(alleles_main), ")")
    val2 <- "-"
    val3 <- paste0(length(genes_unknown), " (", length(alleles_unknown), ")")
    val3[val3 == "0 (0)"] <- "-"
    data.frame(
      "IGHV Family" = chunk_fam$Family[1],
      "F genes (alleles)" = val1,
      "Open reading frame genes (alleles)" = val2,
      "F genes (alleles) unknown" = val3,
      check.names = FALSE,
      stringsAsFactors = FALSE)
  }))
  result
}

make_table1(genbank_alleles)

# groups <- c(
#   "F known" = "Functional main",
#   "F sister" = "Functional sister",
#   "F unknown" = "Functional NA",
#   "NF known" = "Non-functional main",
#   "ORF known" = "ORF main",
#   "ORF unknown" = "ORF NA",
#   "other known" = "NA main",
#   "other sister" = "NA sister",
#   "other unknown" = "NA unknown",
#   "other" = "NA NA")
# genbank_alleles$Group <- factor(
#   with(genbank_alleles, paste(GeneCategory, ScaffoldLocusGroup)),
#   levels = unname(groups),
#   labels = names(groups))
#
# tables <- list()
# tables$S1A <- with(subset(genbank_alleles, ! AccessionDescriptionPartial & Segment == "IGHV"), table(Family, droplevels(Group)))
# tables$S1B <- with(subset(genbank_alleles, ! AccessionDescriptionPartial & Segment == "IGHD"), table(Family))
# tables$S1C <- with(subset(genbank_alleles, ! AccessionDescriptionPartial & Segment == "IGKV"), table(Family, droplevels(Group)))
# 
# lapply(split(genbank_alleles, genbank_alleles$Family), function(chunk) {
#   primary <- subset(chunk, ! is.na(GeneCategory))
#   secondary <- subset(chunk, is.na(GeneCategory))
#   })
