parse_groupings <- function(txt) {
  result <- data.frame(
    Allele = txt,
    stringsAsFactors = FALSE)
  result$Gene <- sub("\\*.*$", "", result$Allele)
  result$Family <- sub("-.*$", "", result$Gene)
  result$Segment <- sub("^(IG[HLK].).*$", "\\1", result$Family)
  result$Locus <- substr(result$Gene, 1, 3)
  result
}

parse_groupings_sonar <- function(txt) {
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
  # Looks like these are all gaps.  Drop them first.
  genbank <- subset(genbank, feature_qualifier_gene != "")

  # There are one or more rows associated with any particular gene, including
  # the genomic sequence, maybe CDS, maybe RSS, etc.  We'll collapse it down to
  # exactly one row per allele.
  genes <- do.call(rbind, lapply(split(genbank, genbank$feature_qualifier_gene), function(chunk) {
    # At this point every chunk should have a row for feature_type "gene" we'll
    # start with
    chunk_out <- subset(
      chunk,
      feature_type == "gene",
      select = c(feature_qualifier_gene, gbf_id, gbf_description, gbf_seqlen, feature_seq))
    colnames(chunk_out) <- c("AlleleOrig", "Accession", "AccessionDescription", "GBFLen", "SeqGenomic")
    # Add separate sequence columns for CDS and translation, if present.
    cds <- subset(chunk, feature_type == "CDS")
    if (nrow(cds) == 0) {
      seq_cds <- ""
      seq_aa <- ""
      product <- ""
    } else if (nrow(cds) == 1) {
      seq_cds <- cds$feature_seq
      seq_aa <- cds$feature_qualifier_translation
      product <- cds$feature_qualifier_product
    } else {
      stop("multiple CDSs?")
    }
    chunk_out$SeqCDS <- seq_cds
    chunk_out$SeqAA <- seq_aa
    chunk_out$Product <- product
    chunk_out
  }))

  # Split out the ontological stuff
  genes <- cbind(genes, parse_groupings(genes$Allele))
  if (all(genes$AlleleOrig == genes$Allele)) {
    genes <- subset(genes, select = -c(AlleleOrig))
  } else {
    stop("unexpected allele name")
  }
  # Also parse out details from the accession descriptions
  genes$AccessionDescriptionPartial <- grepl("partial (cds|sequence)", genes$AccessionDescription)
  genes$AccessionDescriptionFunctional <- NA
  genes$AccessionDescriptionFunctional[grepl("nonfunctional", genes$AccessionDescription)] <- FALSE
  # It looks like we can distinguish the two sequencing approaches by how the
  # GenBank entries are labeled.
  genes$Dataset <- factor(
    genes$AccessionDescriptionPartial,
    levels = c(FALSE, TRUE),
    labels = c("WGS", "Targeted"))
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

load_paper <- function() {
  sheets <- c(paste0("fig", 1:6), paste0("suppsheet", 1:3))
  names(sheets) <- sheets
  paper <- lapply(sheets, function(thing) {
    read.csv(
      file.path("from-paper", paste0(thing, ".csv")),
      stringsAsFactors = FALSE)
  })
  paper
}

parse_genes <- function(paper) {
  paper$fig1$Gene <- paste0("IGHV", paper$fig1$Gene)
  paper$fig2$Gene <- paste0("IGKV", paper$fig2$Gene)
  paper$fig3$Gene <- paste0("IGLV", paper$fig3$Gene)
  part1 <- do.call(
    rbind, lapply(paper[c("fig1", "fig2", "fig3")], function(sheet) {
    if (! "LocusGroup" %in% colnames(sheet)) {
      sheet$LocusGroup <- as.character(NA)
    }
    sheet[, c("Gene", "Category", "LocusGroup")]
  }))
  colnames(part1)[1] <- "GeneOrig"
  part1 <- cbind(part1, parse_groupings(part1$GeneOrig))
  part1 <- part1[
    , -match(c("GeneOrig", "Allele"), colnames(part1))]
  
  part2 <- do.call(
    rbind, lapply(paper[c("fig4", "fig5", "fig6")], function(sheet) {
      sheet$Category <- "Functional"
      if (! "LocusGroup" %in% colnames(sheet)) {
        sheet$LocusGroup <- as.character(NA)
      }
      sheet[, c("Allele", "Category")]
    }))
  colnames(part2)[1] <- "AlleleOrig"
  part2 <- cbind(part2, parse_groupings(part2$AlleleOrig))
  part2 <- part2[
    , -match(c("Allele", "AlleleOrig"), colnames(part2))]
  part2$LocusGroup <- as.character(NA)
  result <- rbind(part1, subset(part2, ! Gene %in% part1$Gene))
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

merge_metadata <- function(genbank_alleles, metadata) {
  idx <- match(genbank_alleles$Gene, metadata$genes$Gene)
  genbank_alleles$GeneCategory <- metadata$genes$Category[idx]
  genbank_alleles$GeneLocusGroup <- metadata$genes$LocusGroup[idx]
  idx <- match(genbank_alleles$GBFLen, metadata$scaffolds$Length)
  genbank_alleles$Scaffold <- metadata$scaffolds$Scaffold[idx]
  genbank_alleles$ScaffoldLocusGroup <- metadata$scaffolds$LocusGroup[idx]
  genbank_alleles
}
