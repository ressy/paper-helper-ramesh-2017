
# Helpers -----------------------------------------------------------------


parse_groupings <- function(txt) {
  result <- data.frame(
    Gene = sub("^ORF_", "", sub("\\*.*$", "", txt)),
    stringsAsFactors = FALSE)
  result$Family <- sub("-.*$", "", result$Gene)
  result$Segment <- sub("^(IG[HLK].).*$", "\\1", result$Family)
  result$Locus <- substr(result$Gene, 1, 3)
  result
}

parse_genbank_genes <- function(genbank) {
  genes <- subset(
    genbank,
    feature_type == "gene",
    select = c(feature_qualifier_gene, gbf_id, gbf_description, feature_seq))
  colnames(genes) <- c("Allele", "Accession", "AccessionDescription", "Seq")
  # Split out the ontological stuff
  genes <- cbind(genes, parse_groupings(genes$Allele))
  # Also parse out details from the accession descriptions
  genes$Partial <- grepl("partial (cds|sequence)", genes$AccessionDescription)
  genes$Functional <- NA
  genes$Functional[grepl("nonfunctional", genes$AccessionDescription)] <- FALSE
  genes
}

load_sonar_alleles <- function() {
  sonar_alleles <- do.call(
    rbind,
    lapply(
      list.files(
        "SONAR/germDB", pattern="Ig(HD|HKL[VDJ])_BU_DD\\.fasta$",
        full.names = TRUE), dnar::read.fa))
  colnames(sonar_alleles) <- c("Allele", "Seq")
  sonar_alleles <- cbind(sonar_alleles, parse_groupings(sonar_alleles$Allele))
  sonar_alleles
}


# Go ----------------------------------------------------------------------


genbank <- read.csv("converted/all.csv", stringsAsFactors = FALSE)
genbank_alleles <- parse_genbank_genes(genbank)
sonar_alleles <- load_sonar_alleles()

# Across the board the SONAR set has fewer entries for each segment than the 
# ones in GenBank from the paper.  Also, some of the names don't match up, with
# an extra suffix on some of the SONAR ones.
segments <- c(paste0("IGH", c("V", "D", "J")), "IGLV", "IGLJ", "IGKV", "IGKJ")
table(subset(sonar_alleles, Segment %in% segments)$Segment)[segments]
table(subset(genbank_alleles, Segment %in% segments & ! Partial)$Segment)[segments]
