
# Helpers -----------------------------------------------------------------


parse_groupings <- function(txt) {
  result <- data.frame(
    Prefix = sub("_?IG.*", "", txt),
    Suffix = sub(
      "(.*IG[HKL]V[0-9]+-[^-]+-?|.*IG.[^V].*)", "", sub("\\*.*", "", txt)),
    Allele = sub("(?:[^_]*_)", "", sub("(.*IG.[V].*)-.\\*", "\\1*", txt)),
    stringsAsFactors = FALSE)
  result$Gene <- sub("\\*.*$", "", result$Allele)
  result$Family <- sub("-.*$", "", result$Gene)
  result$Segment <- sub("^(IG[HLK].).*$", "\\1", result$Family)
  result$Locus <- substr(result$Gene, 1, 3)
  result
}

parse_genbank_genes <- function(genbank) {
  genes <- subset(
    genbank,
    feature_type == "gene",
    select = c(feature_qualifier_gene, gbf_id, gbf_description, gbf_seqlen, feature_seq))
  colnames(genes) <- c("AlleleOrig", "Accession", "AccessionDescription", "GBFLen", "Seq")
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
  colnames(sonar_alleles) <- c("AlleleOrig", "Seq")
  sonar_alleles <- cbind(sonar_alleles, parse_groupings(sonar_alleles$AlleleOrig))
  sonar_alleles
}


# Go ----------------------------------------------------------------------


genbank <- read.csv("converted/all.csv", stringsAsFactors = FALSE)
genbank_alleles <- parse_genbank_genes(genbank)
sonar_alleles <- load_sonar_alleles()

# SONAR generally has fewer entries than what the paper described, but not for
# all (for example IGKV), and it includes ones labeled "partial" that the paper
# seems to exclude from the final counts.
segments <- c(paste0("IGH", c("V", "D", "J")), "IGLV", "IGLJ", "IGKV", "IGKJ")
table(subset(genbank_alleles, Segment %in% segments)$Segment)[segments]
table(subset(genbank_alleles, Segment %in% segments & ! Partial)$Segment)[segments]
table(subset(sonar_alleles, Segment %in% segments)$Segment)[segments]

sheets <- c(paste0("fig", 1:6), paste0("suppsheet", 1:3))
names(sheets) <- sheets
paper <- lapply(sheets, function(thing) {
  read.csv(
    file.path("from-paper", paste0(thing, ".csv")),
    stringsAsFactors = FALSE)
})

parse_genes <- function(paper) {
  # paper$fig1$Locus <- "IGH"
  # paper$fig2$Locus <- "IGK"
  # paper$fig3$Locus <- "IGL"
  paper$fig1$Gene <- paste0("IGHV", paper$fig1$Gene)
  paper$fig2$Gene <- paste0("IGKV", paper$fig2$Gene)
  paper$fig3$Gene <- paste0("IGLV", paper$fig3$Gene)
  result <- do.call(
    rbind, lapply(paper[c("fig1", "fig2", "fig3")], function(sheet) {
    if (! "LocusGroup" %in% colnames(sheet)) {
      sheet$LocusGroup <- as.character(NA)
    }
    sheet[, c("Gene", "Category", "LocusGroup")]
  }))
  colnames(result)[1] <- "GeneOrig"
  result <- cbind(result, parse_groupings(result$GeneOrig))
  result <- result[
    , -match(c("GeneOrig", "Allele", "Prefix", "Suffix"), colnames(result))]
  rownames(result) <- NULL
  result
}

parse_alleles <- function(paper) {
  # paper$fig4$Allele <- paste0("IGHV", paper$fig4$Allele)
  # paper$fig5$Allele <- paste0("IGKV", paper$fig5$Allele)
  # paper$fig6$Allele <- paste0("IGLV", paper$fig6$Allele)
  result <- do.call(
    rbind, lapply(paper[c("fig4", "fig5", "fig6")], function(sheet) {
      sheet$Category <- "Functional"
      if (! "LocusGroup" %in% colnames(sheet)) {
        sheet$LocusGroup <- as.character(NA)
      }
      sheet[, c("Allele", "Category")]
    }))
  colnames(result)[1] <- "AlleleOrig"
  result <- cbind(result, parse_groupings(result$AlleleOrig))
  result <- result[
    , -match(c("AlleleOrig", "Prefix", "Suffix"), colnames(result))]
  rownames(result) <- NULL
  result
}

parse_scaffolds <- function(paper) {
  paper$suppsheet1$Locus <- "IGH"
  paper$suppsheet2$Locus <- "IGK"
  paper$suppsheet3$Locus <- "IGL"
  result <- with(paper, rbind(suppsheet1, suppsheet2, suppsheet3))
  result
}

metadata <- list(
  genes = parse_genes(paper),
  alleles = parse_alleles(paper),
  scaffolds = parse_scaffolds(paper)
)

idx <- match(genbank_alleles$GBFLen, metadata$scaffolds$Length)
genbank_alleles$Scaffold <- metadata$scaffolds$Scaffold[idx]
genbank_alleles$ScaffoldGroup <- metadata$scaffolds$LocusGroup[idx]
idx <- match(genbank_alleles$Gene, with(metadata$genes, paste0(Locus, "V", Gene)))
genbank_alleles$Category <- metadata$genes$Category[idx]

groups <- c(
  "F known" = "Functional main",
  "F sister" = "Functional sister",
  "F unknown" = "Functional NA",
  "NF known" = "Non-functional main",
  "ORF known" = "ORF main",
  "ORF unknown" = "ORF NA",
  "other known" = "NA main",
  "other sister" = "NA sister",
  "other unknown" = "NA unknown",
  "other" = "NA NA")
genbank_alleles$Group <- factor(
  with(genbank_alleles, paste(Category, ScaffoldGroup)),
  levels = unname(groups),
  labels = names(groups))

tables <- list()
tables$S1A <- with(subset(genbank_alleles, ! Partial & Segment == "IGHV"), table(Family, droplevels(Group)))
tables$S1B <- with(subset(genbank_alleles, ! Partial & Segment == "IGHD"), table(Family))
tables$S1C <- with(subset(genbank_alleles, ! Partial & Segment == "IGKV"), table(Family, droplevels(Group)))

