# TODO: include misc_feature notes!! They explain more.  For example, IGHE-P*01.

# Parse antibody attributes from the Product column text.
parse_fields_from_product <- function(txt) {
  # yes some of the GenBank records do have "lamba" in them.  I'm in no position
  # to criticize considering how often I write "lamda".
  match <- regexec(
    "^immunoglobulin (heavy|lambda|lamba|kappa)(?: chain)? ?(alpha|delta|epsilon|gamma|mu)?(.*)(domain|region) ?([0-9]*)$",
    txt)
  fields <- do.call(rbind, lapply(regmatches(txt, match), `[`, 2:6))
  fields[is.na(fields)] <- ""
  # the apply function has no drop=F apparently.
  for (idx in 1:ncol(fields)) {
    fields[, idx] <- gsub(" *$", "", gsub("^ *", "", fields[, idx]))
  }
  fields <- as.data.frame(fields, stringsAsFactors = FALSE)
  colnames(fields) <- c(
    "Locus", "Class", "region_segment", "region_domain", "Subclass")
  fields$Locus <- c(heavy="IGH", kappa="IGK", lambda="IGL")[fields$Locus]
  fields$Locus[is.na(fields$Locus)] <- ""
  fields$Class[fields$Class == "lamba"] <- "lambda"
  fields$Region <- ifelse(
    fields$region_segment %in% c("variable", "diversity", "joining", "constant"),
    fields$region_segment,
    ifelse(fields$region_domain == "domain", "constant", ""))
  fields$Domain <- ifelse(
    fields$region_domain == "domain", fields$region_segment, "")
  fields <- subset(fields, select = -c(region_domain, region_segment))
  colnames(fields) <- paste0("Product", colnames(fields))
  fields
}

# Parse antibody attributes from the AccessionDescription column text.
parse_fields_from_accession_description <- function(txt) {
  match <- regexec(
    "^Macaca mulatta nonfunctional immunoglobulin (heavy|lambda|kappa) chain (alpha|delta|epsilon|gamma|mu)? ?(.*) (domain|region) ?([0-9]*)",
    txt)
  fields <- do.call(rbind, lapply(regmatches(txt, match), `[`, 2:6))
  fields[is.na(fields)] <- ""
  fields <- as.data.frame(fields, stringsAsFactors = FALSE)
  colnames(fields) <- c(
    "Locus", "Class", "region_segment", "region_domain", "Subclass")
  fields$Locus <- c(heavy="IGH", kappa="IGK", lambda="IGL")[fields$Locus]
  fields$Locus[is.na(fields$Locus)] <- ""
  fields$Region <- ifelse(
    fields$region_segment %in% c("variable", "diversity", "joining", "constant"),
    fields$region_segment,
    ifelse(fields$region_domain == "domain", "constant", ""))
  fields$Domain <- ifelse(
    fields$region_domain == "domain", fields$region_segment, "")
  fields <- subset(fields, select = -c(region_domain, region_segment))
  colnames(fields) <- paste0("AccessionDescription", colnames(fields))
  fields
}

# Parse antibody attributes from the Allele column text.
parse_fields_from_allele <- function(txt, prefix="Allele") {
  fields <- data.frame(
    Allele = txt,
    stringsAsFactors = FALSE)
  fields$Gene <- sub("\\*.*$", "", fields$Allele)
  fields$Family <- sub("-.*$", "", fields$Gene)
  fields$Segment <- sub("^(IG[HLK].).*$", "\\1", fields$Family)
  fields$Locus <- substr(fields$Gene, 1, 3)
  class_gene_lut <- c(
    IGHA="alpha",
    IGHD="delta",
    IGHE="epsilon",
    IGHEP="epsilon",
    "IGHE-P"="epsilon",
    IGHG="gamma",
    IGHG1="gamma",
    IGHG2="gamma",
    IGHG3="gamma",
    IGHG4="gamma",
    IGHM="mu")
  subclass_gene_lut <- c(IGHG1=1, IGHG2=2, IGHG3=3, IGHG4=4)
  fields$Class <- class_gene_lut[fields$Gene]
  fields$Class[is.na(fields$Class)] <- ""
  fields$Subclass <- subclass_gene_lut[fields$Gene]
  fields$Subclass[is.na(fields$Subclass)] <- ""
  region_segment_lut <- c(
    IGHV = "variable",
    IGKV = "variable",
    IGLV = "variable",
    IGHD = "diversity", # needs overriding below for one name clash
    IGHJ = "joining",
    IGKJ = "joining",
    IGLJ = "joining",
    IGHA = "constant",
    IGHE = "constant",
    IGHG = "constant",
    IGHM = "constant",
    IGKC = "constant",
    IGLC = "constant")
  fields$Region <- region_segment_lut[fields$Segment]
  fields$Region[fields$Class != ""] <- "constant"
  fields$Domain <- gsub("^[^_]*_?", "", fields$Allele)
  # I think Family and Segment might only be applicable for VDJ.
  vdj <- c("variable", "diversity", "joining")
  fields$Family[! fields$Region %in% vdj] <- ""
  fields$Segment[! fields$Region %in% vdj] <- ""
  colnames(fields) <- paste0(prefix, colnames(fields))
  colnames(fields)[1] <- "Allele"
  fields
}

patch_genbank_entries <- function(genbank) {
  # MF989453.1: kappa versus lambda
  idx <- which(with(genbank, Accession == "MF989453.1" & AlleleOrig == "IGLV2-ABX*01"))
  if (length(idx) == 1) {
    genbank$Product[idx] <- gsub(
      "^immunoglobulin kappa chain variable region$",
      "immunoglobulin lambda chain variable region",
      genbank$Product[idx])
  } else if (length(idx) > 1) {
    stop("More than one IGLV2-ABX*01?")
  }
  # MF989451.1: delta versus mu
  idx <- which(with(genbank, Accession == "MF989451.1" & AlleleOrig == "IGHM*01"))
  if (length(idx) == 1) {
    genbank$Product[idx] <- gsub(
      "^immunoglobulin heavy chain delta constant region$",
      "immunoglobulin heavy chain mu constant region",
      genbank$Product[idx])
  } else if (length(idx) > 1) {
    stop("More than one IGHM*01?")
  }
  genbank
}

parse_genbank_entries <- function(genbank) {
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

  # Apply fixes before trying to parse more details
  genes <- patch_genbank_entries(genes)

  # Split out the ontological stuff
  genes <- cbind(
    genes,
    parse_fields_from_allele(genes$Allele),
    parse_fields_from_accession_description(genes$AccessionDescription),
    parse_fields_from_product(genes$Product))
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
  # Apply an initial ordering
  genes <- genes[with(genes, order(Accession, Allele)), ]
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
  sonar_alleles <- cbind(sonar_alleles, parse_fields_from_allele(sonar_alleles$AlleleOrig))
  sonar_alleles
}

load_paper <- function(dirpath) {
  sheets <- c(paste0("fig", 1:6), paste0("suppsheet", 1:3))
  names(sheets) <- sheets
  paper <- lapply(sheets, function(thing) {
    read.csv(
      file.path(dirpath, paste0(thing, ".csv")),
      stringsAsFactors = FALSE)
  })
  paper
}

parse_paper_genes <- function(paper) {
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
  part1 <- cbind(part1, parse_fields_from_allele(part1$GeneOrig, ""))
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
  part2 <- cbind(part2, parse_fields_from_allele(part2$AlleleOrig, ""))
  part2 <- part2[
    , -match(c("Allele", "AlleleOrig"), colnames(part2))]
  part2$LocusGroup <- as.character(NA)
  result <- rbind(part1, subset(part2, ! Gene %in% part1$Gene))
  rownames(result) <- NULL
  result
}

parse_paper_scaffolds <- function(paper) {
  paper$suppsheet1$Locus <- "IGH"
  paper$suppsheet2$Locus <- "IGK"
  paper$suppsheet3$Locus <- "IGL"
  result <- with(paper, rbind(suppsheet1, suppsheet2, suppsheet3))
  result
}

parse_paper <- function(dirpath) {
  paper <- load_paper(dirpath)
  metadata <- list(
    genes = parse_paper_genes(paper),
    scaffolds = parse_paper_scaffolds(paper)
  )
  metadata
}

# Merge info from the paper with GenBank entries
merge_metadata <- function(genbank_alleles, metadata) {
  idx <- match(genbank_alleles$AlleleGene, metadata$genes$Gene)
  genbank_alleles$GeneCategory <- metadata$genes$Category[idx]
  genbank_alleles$GeneLocusGroup <- metadata$genes$LocusGroup[idx]
  idx <- match(genbank_alleles$GBFLen, metadata$scaffolds$Length)
  genbank_alleles$Scaffold <- metadata$scaffolds$Scaffold[idx]
  genbank_alleles$ScaffoldLocusGroup <- metadata$scaffolds$LocusGroup[idx]
  genbank_alleles
}

# Condense multiple sources of info into one column each.
# A number of things (constant region domain, locus, etc.) can be parsed out
# from a few different places so this collapses them down to one column each.
collapse_fields <- function(genbank_alleles) {
  fields <- c(
    "Gene", "Family", "Segment", "Locus",
    "Class", "Subclass", "Region", "Domain")
  sources <- c("Allele", "Product", "AccessionDescription")
  for (field in fields) {
    columns <- paste0(sources, field)
    columns <- columns[columns %in% colnames(genbank_alleles)]
    genbank_alleles[[field]] <- apply(
      genbank_alleles[, columns, drop=FALSE], 1, function(vec) {
      vec <- unique(vec[vec != "" & ! is.na(vec)])
      if (length(vec) == 0) {
        ""
      } else if (length(vec) == 1) {
        vec
      } else {
        stop(paste0(
          "Mismatch between data sources for ",
          field,
          ": ",
          paste(vec, collapse="/")))
      }
    })
    # remove the original separate columns
    for (column in columns) {
      genbank_alleles[[column]] <- NULL
    }
  }
  genbank_alleles
}

# Order rows and columns consistently
finalize_table <- function(genbank_alleles) {
  columns <- c(
    "Locus",    # IGH, IGK, or IGL
    "Region",   # variable, diversity, joining, or constant
    "Domain",   # for heavy constant region, if the sequence is only one domain
    "Class",    # for heavy constant: alpha, delta, gamma, epsilon, or mu
    "Subclass", # for heavy constant, like 1 for IGHG1
    "Allele",   # Full sequence identifier
    "Gene",     # Everything before the * from Allele
    "Family",   # For VDJ sequences, IG[HKL][VDJ][0-9]P?
    "Segment",  # For VDJ sequences, IG[HKL][VDJ]
    "Accession",
    "AccessionDescription",
    "Product",
    "AccessionDescriptionPartial",
    "AccessionDescriptionFunctional",
    "GBFLen",             # length of genomic seq for each full GenBank entry (across all regions)
    "Dataset",            # WGS or Targeted
    "GeneCategory",       # Functional, Non-functional, or ORF
    "GeneLocusGroup",     # main or sister
    "Scaffold",           # scaffold ID and size
    "ScaffoldLocusGroup", # main, sister, or unknown
    "SeqGenomic",         # Full genomic sequence with introns
    "SeqCDS",             # Coding sequence
    "SeqAA")              # Translation of coding sequence
  if (! all(sort(columns) == sort(colnames(genbank_alleles)))) {
    stop("unexpected columns")
  }
  # Order columns, and then order rows according to columns
  genbank_alleles <- genbank_alleles[, columns]
  genbank_alleles <- genbank_alleles[do.call(order, genbank_alleles), ]
  genbank_alleles
}