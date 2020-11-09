genbank <- read.csv("converted/all.csv", stringsAsFactors = FALSE)

parse_genbank_genes <- function(genbank) {
  genes <- subset(
    genbank,
    feature_type == "gene",
    select = c(feature_qualifier_gene, gbf_id, gbf_description, feature_seq))
  colnames(genes) <- c("Allele", "Accession", "AccessionDescription", "Seq")
  # Split out the ontological stuff
  genes$Gene <- sub("\\*.*$", "", genes$Allele)
  genes$Family <- sub("-.*$", "", genes$Gene)
  genes$Segment <- sub("^(IG[HLK].).*$", "\\1", genes$Family)
  genes$Locus <- substr(genes$Gene, 1, 3)
  # Also parse out details from the accession descriptions
  genes$Partial <- grepl("partial (cds|sequence)", genes$AccessionDescription)
  genes$Functional <- NA
  genes$Functional[grepl("nonfunctional", genes$AccessionDescription)] <- FALSE
  genes
}

genes <- parse_genbank_genes(genbank)