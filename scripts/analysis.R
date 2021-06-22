# Main --------------------------------------------------------------------


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

    # Then parse all the other rows.
    extras <- subset(chunk, ! feature_type %in% c("gene"))
    bits <- with(extras, do.call(paste0, list(
      0 + (feature_type == "CDS"),
      0 + (feature_type == "misc_recomb"),
      0 + (feature_type == "misc_feature"),
      0 + (feature_qualifier_note == "recombination signal sequence"),
      0 + (feature_qualifier_note == "mutated recombination signal sequence"),
      0 + (grepl("^nonfunctional immunoglobulin", feature_qualifier_note)))))
    extras$feature_category <- factor(
      bits,
      levels = c("010100", "001010", "001001", "100000"),
      labels = c("RSS", "RSS", "CDS", "CDS"))
    extras$feature_mutated <- c("010100"=FALSE, "001010"=TRUE, "001001"=TRUE, "100000"=FALSE)[bits]
    
    if (any(is.na(extras$feature_category))) {
      stop(paste("unrecognized GBF features for", chunk_out$AlleleOrig))
    }
    
    if (sum(extras$feature_category == "RSS") > 2) {
      stop(paste("> 2 RSS entries for", chunk_out$AlleleOrig))
    }
    if (sum(extras$feature_category == "CDS") > 1) {
      stop(paste("> 1 CDS entries for", chunk_out$AlleleOrig))
    }
    if (nrow(extras) > 0) {
      # assume order based on row order, but flip first if we're on the other
      # strand
      strand = unique(chunk$feature_strand)
      if (length(strand) > 1) {
        stop(paste("inconsistent strand for", chunk_out$AlleleOrig))
      }
      if (strand == -1) {
        extras <- extras[seq(nrow(extras), 1, -1), ]
      }
    }
      
    extras$feature_category_final <- as.character(extras$feature_category)
    idx <- 1:nrow(extras)
    idx_cds <- match("CDS", extras$feature_category)
    extras$feature_category_final[extras$feature_category == "RSS" & idx < idx_cds] <- "RSSUpstream"
    extras$feature_category_final[extras$feature_category == "RSS" & idx > idx_cds] <- "RSSDownstream"
    extras$feature_category_final <- factor(
      extras$feature_category_final,
      levels = c("RSSUpstream", "CDS", "RSSDownstream"))

    # Add separate sequence columns for CDS and translation, if present.
    cds <- subset(extras, feature_category_final == "CDS")
    if (nrow(cds) == 0) {
      seq_cds <- ""
      seq_aa <- ""
      product <- ""
      cds_mut <- NA
    } else {
      seq_cds <- cds$feature_seq
      seq_aa <- cds$feature_qualifier_translation
      product <- cds$feature_qualifier_product
      cds_mut <- cds$feature_mutated
    }
    rss1 <- subset(extras, feature_category_final == "RSSUpstream")
    if (nrow(rss1) == 0) {
      seq_rss1 <- ""
      rss1_mut <- NA
    } else {
      seq_rss1 <- rss1$feature_seq
      rss1_mut <- rss1$feature_mutated
    }
    rss2 <- subset(extras, feature_category_final == "RSSDownstream")
    if (nrow(rss2) == 0) {
      seq_rss2 <- ""
      rss2_mut <- NA
    } else {
      seq_rss2 <- rss2$feature_seq
      rss2_mut <- rss2$feature_mutated
    }
    
    misc_notes <- sort(subset(chunk, feature_type == "misc_feature")$feature_qualifier_note)
    misc_notes <- paste(misc_notes, collapse="; ")
    chunk_out$MiscNotes <- misc_notes
    chunk_out$RSSUpstream <- seq_rss1
    chunk_out$SeqCDS <- seq_cds
    chunk_out$RSSDownstream <- seq_rss2
    chunk_out$SeqAA <- seq_aa
    chunk_out$Product <- product
    chunk_out$RSSUpstreamMutated <- rss1_mut
    chunk_out$RSSDownstreamMutated <- rss2_mut
    chunk_out$CDSMutated <- cds_mut
    chunk_out
  }))
  
  # Apply fixes before trying to parse more details
  genes <- patch_genbank_entries(genes)
  
  # Split out the ontological stuff
  genes <- cbind(
    genes,
    parse_fields_from_allele(genes$Allele),
    parse_fields_from_accession_description(genes$AccessionDescription),
    parse_fields_from_product(genes$Product),
    parse_fields_from_misc_notes(genes$MiscNotes))
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

# Order rows and columns consistently and make some last adjustments
finalize_table <- function(genbank_alleles) {
  
  # nearly all other pseudogenes don't get a hyphen
  genbank_alleles$Allele[genbank_alleles$Allele == "IGKJ5-P*01"] <- "IGKJ5P*01"
  genbank_alleles$Gene[genbank_alleles$Gene == "IGKJ5-P"] <- "IGKJ5P"
  # also, elsewhere this same gene is referred to as IGHEP, not IGHE-P.
  genbank_alleles$Allele[genbank_alleles$Allele == "IGHE-P*01"] <- "IGHEP*01"
  genbank_alleles$Gene[genbank_alleles$Gene == "IGHE-P"] <- "IGHEP"
  
  # Set up a per-every-allele category, defaulting to the associated
  # GeneCategory if any.
  genbank_alleles$Category <- genbank_alleles$GeneCategory
  
  # IGHV
  # There are four ORF genes with five alleles in table 1 and five uncategorized
  # IGHV entries across four genes with coding reigons defined, all in IGHV3.
  genbank_alleles$Category[
    with(genbank_alleles, Segment == "IGHV" & is.na(Category))] <- "ORF"

  # IGKV
  # There is exactly one ORF gene+allele from the Targeted dataset for IGKV1 in
  # table 3 and one uncategorized IGKV1 left, and ditto for IGKV2, but with 2
  # genes and 5 alleles.  That accounts for all uncategorized IGKV rows.
  genbank_alleles$Category[
    with(genbank_alleles, Segment == "IGKV" & is.na(Category))] <- "ORF"
  
  # IGLV
  # The two IGLV7 ORF entries (table 4) are already labeled which just leaves
  # the three IGLV8 for IGLV8-AEE.  That accounts for all uncategorized IGLV
  # rows.
  genbank_alleles$Category[
    with(genbank_alleles, Segment == "IGLV" & is.na(Category))] <- "ORF"
  
  # Otherwise, label Non-functional if we can infer that from GenBank metadata,
  # or Functional otherwise.
  genbank_alleles$Category <- with(
    genbank_alleles,
    ifelse(
      is.na(Category),
      ifelse(
        Functional == "F" | SeqAA == "" | RSSUpstreamMutated == "T" | RSSDownstreamMutated == "T",
        "Non-functional",
        "Functional"),
      Category))
  
  # In IGHV we have two functional rows that I couldn't match to named scaffolds
  # in the supplemental tables, for IGHV3-AFA and IGHV4-AGP.  If one of those is
  # on the sister locus and one on an unknown locus then the checks should match
  # up.  If IGHV4's is unknown that would match table 1 best which would imply
  # sister locus for the IGHV3 one.
  genbank_alleles$ScaffoldLocusGroup[with(genbank_alleles, Gene == "IGHV4-AGP" & is.na(ScaffoldLocusGroup))] <- "unknown"
  genbank_alleles$ScaffoldLocusGroup[with(genbank_alleles, Gene == "IGHV3-AFA" & is.na(ScaffoldLocusGroup))] <- "sister"
  # With that the only mismatch for IGHV functional genes in table 1 is the
  # "1 (3)" in my IGHV3 row versus the paper's IGHV2 row, but I think that was a
  # typo as it also clashes with figure 4.
  
  columns <- c(
    "Locus",    # IGH, IGK, or IGL
    "Region",   # variable, diversity, joining, or constant
    "Domain",   # for heavy constant region, if the sequence is only one domain
    "Class",    # for heavy constant: alpha, delta, gamma, epsilon, or mu
    "Subclass", # for heavy constant, like 1 for IGHG1
    "Allele",   # Full sequence identifier
    "Gene",     # Everything before the * from Allele
    "Family",   # For VDJ sequences, IG[HKL][VDJ][0-9]
    "Segment",  # For VDJ sequences, IG[HKL][VDJ]
    "Accession",
    "AccessionDescription",
    "Product",
    "MiscNotes",
    "Partial",
    "Functional",
    "GBFLen",               # length of genomic seq for each full GenBank entry (across all regions)
    "Dataset",              # WGS or Targeted
    "GeneCategory",         # Functional, Non-functional, or ORF
    "Category",             # Inferred per-allele category
    "GeneLocusGroup",       # main or sister
    "Scaffold",             # scaffold ID and size
    "ScaffoldLocusGroup",   # main, sister, or unknown
    "SeqGenomic",           # Full genomic sequence with introns
    "RSSUpstream",          # RSS at 5' end, if relevant and if detected
    "SeqCDS",               # Coding sequence
    "RSSDownstream",        # RSS at 3' end, if relevant and if detected
    "RSSUpstreamMutated",   # Is RSSUpstream mutated?
    "CDSMutated",           # Is SeqCDS mutated?
    "RSSDownstreamMutated", # Is RSSDownstream mutated
    "SeqAA")                # Translation of coding sequence
  if (paste(sort(columns), collapse="/") != paste(sort(colnames(genbank_alleles)), collapse="/")) {
    stop(paste0(
      "unexpected columns.\nexpected: ",
      paste(sort(columns), collapse="/"),
      "\nobserved: ",
      paste(sort(colnames(genbank_alleles)), collapse="/")))
  }
  # Order columns, and then order rows according to columns
  genbank_alleles <- genbank_alleles[, columns]
  genbank_alleles <- genbank_alleles[do.call(order, genbank_alleles), ]
  genbank_alleles
}


# Parsing and collapsing fields -------------------------------------------


# Parse antibody attributes from the Product column text.
parse_fields_from_product <- function(txt) {
  fields <- parse_ig_desc(txt)
  colnames(fields) <- paste0("Product", colnames(fields))
  fields
}

# Parse antibody attributes from the AccessionDescription column text.
parse_fields_from_accession_description <- function(txt) {
  fields <- parse_ig_desc(txt)
  colnames(fields) <- paste0("AccessionDescription", colnames(fields))
  fields
}

parse_fields_from_misc_notes <- function(txt) {
  fields <- parse_ig_desc(txt)
  fields$Functional <- ! grepl("nonfunctional.*due to mutation", txt)
  fields$Functional[fields$Functional] <- NA
  colnames(fields) <- paste0("MiscNotes", colnames(fields))
  fields
}

# Parse antibody attributes from the Allele column text.
parse_fields_from_allele <- function(txt, prefix="Allele") {
  fields <- data.frame(
    Allele = txt,
    stringsAsFactors = FALSE)
  fields$Gene <- sub("\\*.*$", "", fields$Allele)
  fields$Family <- sub("P$", "", sub("-.*$", "", fields$Gene))
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

# There are a few different places where we run into freeform descriptions like:
# "immunoglobulin lambda chain joining region"
# "nonfunctional immunoglobulin heavy chain epsilon constant region due to mutation"
# etc.
# This will return a data frame of the common fields.
parse_ig_desc <- function(txt) {
  # yes some of the GenBank records do have "lamba" in them.  I'm in no position
  # to criticize considering how often I write "lamda".
  match <- regexec(
    "immunoglobulin (heavy|lambda|lamba|kappa)(?: chain)? ?(alpha|delta|epsilon|gamma|mu)? ?(constant|variable|diversity|joining|[-CHS0-9]+) ?(domain|region) ?([0-9]*)",
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
  fields
}

# Condense multiple sources of info into one column each.
# A number of things (constant region domain, locus, etc.) can be parsed out
# from a few different places so this collapses them down to one column each and
# removes the separate ones.
collapse_fields <- function(genbank_alleles) {
  fields <- c(
    "Gene", "Family", "Segment", "Locus",
    "Class", "Subclass", "Region", "Domain")
  sources <- c("Allele", "Product", "AccessionDescription", "MiscNotes")
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
  
  # logicals
  fields <- c("Functional", "Partial",
              "RSSUpstreamMutated", "CDSMutated", "RSSDownstreamMutated")
  sources <- c("AccessionDescription", "MiscNotes", "")
  for (field in fields) {
    columns <- paste0(sources, field)
    columns <- columns[columns %in% colnames(genbank_alleles)]
    genbank_alleles[[field]] <- apply(
      genbank_alleles[, columns, drop=FALSE], 1, function(vec) {
        vec <- unique(vec[vec != "" & ! is.na(vec)])
        if (length(vec) == 0) {
          ""
        } else if (length(vec) == 1 && vec == TRUE) {
          "T"
        } else if (length(vec) == 1 && vec == FALSE) {
          "F"
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
      if (! column %in% fields) {
        genbank_alleles[[column]] <- NULL
      }
    }
  }
  
  genbank_alleles
}


# Parsing info from paper -------------------------------------------------


load_paper <- function(dirpath) {
  sheets <- c(
    paste0("fig", 1:6),
    paste0("suppsheet", 1:3),
    paste0("table", c(1, 3, 4)))
  names(sheets) <- sheets
  paper <- lapply(sheets, function(thing) {
    read.csv(
      file.path(dirpath, paste0(thing, ".csv")),
      check.names = FALSE,
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

parse_paper <- function(paper) {
  metadata <- list(
    genes = parse_paper_genes(paper),
    scaffolds = parse_paper_scaffolds(paper)
  )
  metadata
}


# Checks ------------------------------------------------------------------


# Use the final data frame with all info to build a list of tables to sanity
# check everything.
make_check_tables <- function(genbank_alleles) {
  table1 <- make_v_table(subset(genbank_alleles, Segment == "IGHV"))
  table3 <- make_v_table(subset(genbank_alleles, Segment == "IGKV"))
  table4 <- make_v_table(subset(genbank_alleles, Segment == "IGLV"))
  checks <- list(
    table1 = make_final_table1(table1),
    table3 = make_final_table3(table3),
    table4 = make_final_table4(table4),
    checks = make_checks_comparison(genbank_alleles))
  checks
}

# Make one of the IGHV/IGLV/IGKV tables
make_v_table <- function(ramesh_for_segment) {
  output <- do.call(rbind, lapply(unique(ramesh_for_segment$Family), function(family) {
    ramesh_family <- subset(ramesh_for_segment, Family == family)
    orig_genes <- unique(subset(ramesh_family, Dataset == "WGS" & Category == "Functional" & ! ScaffoldLocusGroup %in% "unknown")$Gene)
    orig_genes_unknown <- unique(subset(ramesh_family, Dataset == "WGS" & Category == "Functional" & ScaffoldLocusGroup %in% "unknown")$Gene)
    orig_orf <- unique(subset(ramesh_family, Dataset == "WGS" & Category == "ORF" & SeqAA != "")$Gene)
    extra_genes <- unique(subset(ramesh_family, Dataset == "Targeted" & Category == "Functional" & ! Gene %in% orig_genes & ! Gene %in% orig_genes_unknown)$Gene)
    extra_orf <- unique(subset(ramesh_family, Dataset == "Targeted" & Category == "ORF" & SeqAA != "" & ! Gene %in% orig_orf)$Gene)
    full_alleles <- subset(ramesh_family, Gene %in% c(orig_genes, extra_genes))$Allele
    orf_alleles <- subset(ramesh_family, Gene %in% c(orig_orf, extra_orf))$Allele
    unknown_alleles <- subset(ramesh_family, Gene %in% orig_genes_unknown)$Allele
    data.frame(
      Family = family,
      FunctionalGenesWGS = length(orig_genes),
      FunctionalGenesTargeted = length(extra_genes),
      FunctionalAlleles = length(full_alleles),
      ORFGenesWGS = length(orig_orf),
      ORFGenesTargeted = length(extra_orf),
      ORFAlleles = length(orf_alleles),
      FunctionalGenesWGSUnknown = length(orig_genes_unknown),
      FunctionalAllelesUnknown = length(unknown_alleles),
      stringsAsFactors = FALSE)
  }))
  output <- output[order(as.integer(gsub("[^0-9]", "", output$Family))), ]
  output
}

# just a helper to format text like the paper's tables.
fmt_table_cols <- function(col1, col2, col3) {
  sapply(1:length(col1), function(idx) {
    if (! is.na(col2[idx]) && col2[idx] > 0) {
      txt <- paste(col1[idx], col2[idx], sep = " + ")
    } else {
      txt <- col1[idx]
    }
    txt <- paste0(txt, " (", col3[idx], ")")
    if (txt == "0 (0)") {
      txt <- "-"
    }
    txt
  })
}

make_final_table1 <- function(table1) {
  with(table1, {
    data.frame(
      `IGHV family` = Family,
      `F genes (alleles)` = fmt_table_cols(
        FunctionalGenesWGS, FunctionalGenesTargeted, FunctionalAlleles),
      `Open reading frame genes (alleles)` = fmt_table_cols(
        ORFGenesWGS, ORFGenesTargeted, ORFAlleles),
      `F genes (alleles)` = fmt_table_cols(
        FunctionalGenesWGSUnknown, 0, FunctionalAllelesUnknown),
      stringsAsFactors = FALSE, check.names = FALSE)
  })
}

make_final_table3 <- function(table3) {
  with(table3, {
    data.frame(
      `IGKV family` = Family,
      `F genes (alleles)` = fmt_table_cols(
        FunctionalGenesWGS, FunctionalGenesTargeted, FunctionalAlleles),
      `ORF genes (alleles)` = fmt_table_cols(
        ORFGenesWGS, ORFGenesTargeted, ORFAlleles),
      stringsAsFactors = FALSE, check.names = FALSE)
  })
}

make_final_table4 <- function(table4) {
  with(table4, {
    data.frame(
      `IGLV family` = Family,
      `F genes (alleles)` = fmt_table_cols(
        FunctionalGenesWGS, FunctionalGenesTargeted, FunctionalAlleles),
      `ORF genes (alleles)` = fmt_table_cols(
        ORFGenesWGS, ORFGenesTargeted, ORFAlleles),
      stringsAsFactors = FALSE, check.names = FALSE)
  })
}

# Compare what I have in my final output with various things from the text of
# the paper.
make_checks_comparison <- function(genbank_alleles) {
  options(stringsAsFactors = FALSE)
  wgs <- subset(genbank_alleles, Dataset == "WGS")
  targeted <- subset(genbank_alleles, Dataset == "Targeted")
  
  # a special case for two rows
  known_loc_ighv_genes <- with(wgs, Gene[Segment == "IGHV" & Category == "Functional" & ! ScaffoldLocusGroup %in% "unknown"])
  extra_ighv_genes <- unique(with(targeted, Gene[Segment == "IGHV" & Category == "Functional" & ! Gene %in% wgs$Gene]))
  known_loc_ighv_genes <- c(known_loc_ighv_genes, extra_ighv_genes)
  known_loc_ighv_alleles <- with(genbank_alleles, Allele[Gene %in% known_loc_ighv_genes])
  
  wgs_checks <- with(wgs, do.call(rbind, lapply(list(
    # "IGH: 29 scaffolds were positive for IGH genes. The largest of these is
    # 1.4 Mb long and contains 55 IGHV genes and all 39 IGHD, all 9 IGHJ, and all
    # 9 IGHC genes.
    data.frame(
      Category = "Largest scaffold IGHV genes",
      Expected = 55,
      Observed = sum(Scaffold %in% "scaffold_297" & Segment %in% "IGHV")),
    data.frame(
      Category = "Largest scaffold IGHD genes",
      Expected = 39,
      Observed = sum(Scaffold %in% "scaffold_297" & Segment %in% "IGHD")),
    data.frame(
      Category = "Largest scaffold IGHJ genes",
      Expected = 9,
      Observed = sum(Scaffold %in% "scaffold_297" & Segment %in% "IGHJ")),
    data.frame(
      Category = "Largest scaffold IGHC genes",
      Expected = 9,
      Observed = sum(Scaffold == "scaffold_297" & Locus == "IGH" & Region == "constant"),
      Comments = '"Only 8 functional IGHC genes defined elsewhere; ""Eight IGHC genes were identified including one IGHA gene and four IGHG genes. An additional pseudogene (IGHEP) was found on an unplaced scaffold."" and ""Macaques appear to have only eight IGHC genes"""'),
    # A total of 178 IGHV genes (71 functional) belonging to 7 IGHV families
    # were identified across all scaffolds."
    data.frame(
      Category = "Total IGHV genes",
      Expected = 178,
      Observed = sum(Segment == "IGHV")),
    data.frame(
      Category = "Functional IGHV genes",
      Expected = 71,
      Observed = sum(Segment == "IGHV" & Category %in% "Functional")),
    # IGK: Eight scaffolds were positive for IGK genes. The largest of these is
    # 744 kb long and contains 61 IGKV genes. A total of 119 IGKV (65 functional),
    # 5 IGKJ genes and 1 IGKC gene were identified.
    data.frame(
      Category = "Scaffolds containing IGK genes",
      Expected = 8,
      Observed = length(unique(subset(wgs, Locus == "IGK")$Scaffold))),
    data.frame(
      Category = "Largest scaffold IGKV genes",
      Expected = 61,
      Observed = sum(Scaffold %in% "scaffold1|size744068" & Segment == "IGKV")),
    data.frame(
      Category = "Total IGKV genes",
      Expected = 119,
      Observed = nrow(subset(wgs, Segment == "IGKV"))),
    data.frame(
      Category = "Functional IGKV genes",
      Expected = 65,
      Observed = nrow(subset(wgs, Segment == "IGKV" & Category == "Functional"))),
    data.frame(
      Category = "Total IGKJ genes",
      Expected = 5,
      Observed = nrow(subset(wgs, Locus == "IGK" & Segment == "IGKJ"))),
    data.frame(
      Category = "Total IGKC genes",
      Expected = 1,
      Observed = nrow(subset(wgs, Locus == "IGK" & Region == "constant"))),
    # "IGL: Three scaffolds were positive for IGL genes. The largest of these is
    # 800 kb long and contains 54 IGLV genes, and all seven pairs of alternating
    # IGLJ and IGLC genes."
    data.frame(
      Category = "Scaffolds containing IGL genes",
      Expected = 3,
      Observed = length(unique(subset(wgs, Locus == "IGL")$Scaffold))),
    data.frame(
      Category = "Largest scaffold IGLV genes",
      Expected = 54,
      Observed =  sum(Scaffold %in% "scaffold1|size815129" & Segment == "IGLV")),
    data.frame(
      Category = "Largest scaffold IGLJ genes",
      Expected = 7,
      Observed = sum(Scaffold %in% "scaffold1|size815129" & Segment == "IGLJ")),
    data.frame(
      Category = "Largest scaffold IGLC genes",
      Expected = 7,
      Observed = sum(Scaffold %in% "scaffold1|size815129" & Locus == "IGL" & Region == "constant")),
    data.frame(
      Category = "Total IGLJ genes",
      Expected = 7,
      Observed = nrow(subset(wgs, Locus == "IGL" & Segment == "IGLJ"))),
    data.frame(
      Category = "Total IGLC genes",
      Expected = 7,
      Observed = nrow(subset(wgs, Locus == "IGL" & Region == "constant"))),
    # "A total of 105 IGLV (46 functional, 1 ORF, 58 NF) were identified across
    # scaffolds."
    data.frame(
      Category = "Total IGLV genes",
      Expected = 105,
      Observed = nrow(subset(wgs, Segment == "IGLV"))),
    data.frame(
      Category = "Functional IGLV genes",
      Expected = 46,
      Observed = nrow(subset(wgs, Segment == "IGLV" & Category == "Functional"))),
    data.frame(
      Category = "ORF IGLV genes",
      Expected = 1,
      Observed = nrow(subset(wgs, Segment == "IGLV" & Category == "ORF"))),
    data.frame(
      Category = "Non-functional IGLV genes",
      Expected = 58,
      Observed = nrow(subset(wgs, Segment == "IGLV" & Category == "Non-functional"))),
    # "IGHV: After ordering the scaffolds, we found that 44 functional IGHV genes
    # localized to the main locus, 20 functional genes to the sister locus, and 7
    # functional genes to scaffolds that could not be placed (Figure S1A and Table
    # S1A in Supplementary Material)."
    data.frame(
      Category = "Main locus functional IGHV genes",
      Expected = 44,
      Observed = sum(ScaffoldLocusGroup %in% "main" & Segment %in% "IGHV" & Category %in% "Functional")),
    data.frame(
      Category = "Sister locus functional IGHV genes",
      Expected = 20,
      Observed = sum(ScaffoldLocusGroup %in% "sister" & Segment %in% "IGHV" & Category %in% "Functional")),
    data.frame(
      Category = "Unknown locus functional IGHV genes",
      Expected = 7,
      Observed = sum(ScaffoldLocusGroup %in% "unknown" & Segment %in% "IGHV" & Category %in% "Functional")),
    # "IGHD, J, and C: All the IGHD, IGHJ, and IGHC genes were assembled onto one
    # scaffold, along with 55 functional IGHV genes. Thirty-nine functional IGHD
    # genes belonging to six families were identified— of these, ten formed five
    # pairs of identical genes, and three formed one triplet of identical genes
    # (Table S1B in Supplementary Material). These genes are organized into seven
    # homologous clusters, a cluster prototypically containing one gene from each
    # of the IGHD1-6 families. The observed clusters differ from the prototype by
    # deletions and duplications. Nine IGHJ genes (six functional) were
    # identified, with IGHJ5 duplicated. Eight IGHC genes were identified
    # including one IGHA gene and four IGHG genes. An additional pseudogene
    # (IGHEP) was found on an unplaced scaffold."
    data.frame(
      Category = "IGHD/IGHJ/IGHC genes not on largest scaffold",
      Expected = 0,
      Observed = sum(
        Locus == "IGH" &
          Region %in% c("diversity", "joining", "constant") &
          Gene != "IGHEP" &
          ! Scaffold %in% "scaffold_297")),
    data.frame(
      Category = "Functional IGHD genes",
      Expected = 39,
      Observed = sum(Segment == "IGHD" & Category %in% "Functional")),
    data.frame(
      Category = "Total IGHD families",
      Expected = 6,
      Observed = length(unique(subset(wgs, Segment == "IGHD")$Family))),
    data.frame(
      Category = "Total IGHJ genes",
      Expected = 9,
      Observed = sum(Segment == "IGHJ")),
    data.frame(
      Category = "Functional IGHJ genes",
      Expected = 6,
      Observed = sum(Segment == "IGHJ" & Category %in% "Functional")),
    data.frame(
      Category = "Total IGHC genes",
      Expected = 8,
      # leaving out the psueogene and the separate per-domain entries
      Observed = sum(Locus == "IGH" & Region == "constant" & Domain == "" & Gene != "IGHEP")),
    data.frame(
      Category = "Total IGHA genes",
      Expected = 1,
      Observed = sum(Locus == "IGH" & Region == "constant" & Domain == "" & grepl("IGHA", Gene))),
    data.frame(
      Category = "Total IGHG genes",
      Expected = 4,
      Observed = sum(Locus == "IGH" & Region == "constant" & Domain == "" & grepl("IGHG", Gene))),
    data.frame(
      Category = "Unplaced scaffold IGHEP genes",
      Expected = 1,
      Observed = sum(Gene == "IGHEP" & ScaffoldLocusGroup %in% "unknown")),
    # "IGKV: A total of 119 IGKV genes (Table S2 in Supplementary Material)
    # belonging to 6 IGKV families were identified on four scaffolds constituting
    # the main IGK locus and four contigs constituting the sister locus. We
    # identified 54 functional IGKV genes belonging to six families (Figure S1B in
    # Supplementary Material) and 11 functional genes on the sister IGK locus
    # (Table S2 in Supplementary Material).
    # (119 IGKV already handled above)
    data.frame(
      Category = "IGKV families",
      Expected = 6,
      Observed = length(unique(subset(wgs, Segment == "IGKV")$Family))),
    data.frame(
      Category = "Main locus functional IGKV functional genes",
      Expected = 54,
      Observed = sum(Segment == "IGKV" & Category == "Functional" & ScaffoldLocusGroup %in% "main")),
    data.frame(
      Category = "Main locus IGKV families",
      Expected = 6,
      Observed = length(unique(subset(wgs, Segment == "IGKV" & ScaffoldLocusGroup %in% "main")$Family))),
    data.frame(
      Category = "Sister locus functional IGKV genes",
      Expected = 11,
      Observed = sum(ScaffoldLocusGroup %in% "sister" & Segment == "IGKV" & Category %in% "Functional")),
    # "IGKJ, IGKC: The IGKJ cluster containing four functional genes and one NF
    # gene, and the single IGKC were assembled onto one scaffold, along with 10
    # IGKV genes."
    data.frame(
      Category = "IGKC-containing scaffold functional IGKJ genes",
      Expected = 4,
      Observed = sum(Locus == "IGK" & Segment == "IGKJ" & Scaffold %in% "scaffold5|size145413" & Category %in% "Functional")),
    data.frame(
      Category = "IGKC-containing scaffold non-functional IGKJ genes",
      Expected = 1,
      Observed = sum(Segment == "IGKJ" & Scaffold %in% "scaffold5|size145413" & Category %in% "Non-functional")),
    data.frame(
      Category = "IGKC-containing scaffold IGKV genes",
      Expected = 10,
      Observed = sum(Segment == "IGKV" & Scaffold %in% "scaffold5|size145413")),
    # "IGLV: 105 IGLV genes (Table S3 in Supplementary Material) belonging to 11
    # IGLV families were identified including 46 functional IGLV genes belonging
    # to 10 families (Table S3 and Figure S1C in Supplementary Material) and one
    # ORF belonging to the IGLV7 family. No members of the IGLV9 family were
    # functional."
    # (105 IGLV genes and 11 families handled above)
    data.frame(
      Category = "IGLV functional families",
      Expected = 10,
      Observed = length(unique(subset(wgs, Segment == "IGLV" & Category %in% "Functional")$Family))),
    data.frame(
      Category = "IGLV7 ORF genes",
      Expected = 1,
      Observed = sum(Family == "IGLV7" & Category == "ORF")),
    data.frame(
      Category = "IGLV9 family functional genes",
      Expected = 0,
      Observed = sum(Family == "IGLV9" & Category %in% "Functional")),
    # "IGL-JC: Along with 54 IGLV genes, the IGL-JC cluster was assembled onto one
    # scaffold. Seven pairs of tandem IGLV and IGLC genes were identified. Five of
    # the JC genes were functional. Two of the IGLC genes were identical but were
    # paired with different IGLJ."
    data.frame(
      Category = "Scaffolds containing IGLJ and IGLC genes",
      Expected = 1,
      Observed = length(unique(subset(wgs, Locus == "IGL" & Region %in% c("constant", "joining"))$Scaffold))),
    data.frame(
      Category = "Functional IGLJ genes",
      Expected = 5,
      Observed = sum(Segment == "IGLJ" & Category %in% "Functional")),
    data.frame(
      Category = "Functional IGLC genes",
      Expected = 5,
      Observed = sum(Locus == "IGL" & Region == "constant" & Category %in% "Functional"))),
    function(df) {
      if (! "Comments" %in% colnames(df)) {
        df$Comments <- ""
      }
      df
    })
  ))
  wgs_checks <- cbind(Dataset = "WGS", wgs_checks)
  other_checks <- with(targeted, do.call(rbind, lapply(list(
    # This text matches table 1 but both the text and table disagree with figure
    # 4:
    # "IGHV: 276 Functional/ORF IGHV genes were found in M1–M9, including 142
    # unique functional genes and 5 unique ORFs (Table S5A in Supplementary
    # Material). The largest number of genes/alleles was found in the IGHV3
    # family (Figure 4). The final IGHV allele library, comprising genes from
    # the reference and supplementary assemblies, contains 72 functional genes
    # with 186 allelic forms belonging to 7 families (Table 1), with the IGHV3
    # family being the largest (38 IGHV genes, 112 alleles). In addition, 12
    # alleles were found for 7 functional genes that were identified on contigs
    # that could not be placed (Table 1). 4 IGHV3 ORFs (Table 1) were identified
    # in M1–M9, all of which were absent from the reference IGH locus.
    # Similarity of the genes within each family was highly variable with
    # maximum variation in IGHV3 family."
    # I can't check the 276 or 142 numbers since those overlap between datasets
    # and the full detail of which genes were found in which animals was not
    # given.
    # This one below requires clarification of some of the ScaffoldLocusGroup 
    # entries before it can be correct.
    data.frame(
      Dataset = "Both",
      Category = "Known-location functional IGHV genes",
      Expected = 72,
      Observed = length(known_loc_ighv_genes)),
    data.frame(
      Dataset = "Both",
      Category = "Known-location functional IGHV alleles",
      Expected = 186,
      Observed = length(known_loc_ighv_alleles)),
    data.frame(
      Dataset = "Both",
      Category = "Known-location functional IGHV families",
      Expected = 7,
      Observed = length(unique(subset(genbank_alleles, Segment == "IGHV" & Category == "Functional" & ! ScaffoldLocusGroup %in% "unknown")$Family))),
    # "IGHD: 198 Functional IGHD genes were found in the nine macaques,
    # including 47 unique functional IGHD genes. The IGHD allele library
    # contains the 39 functional genes found in the reference assembly plus 10
    # alleles from the supplementary assemblies."
    # (Same situation as above for the 198 and 47; I can't check those.)
    data.frame(
      Dataset = "WGS",
      Category = "IGHD functional genes",
      Expected = 39,
      Observed = nrow(subset(wgs, Segment == "IGHD" & Category == "Functional"))),
    data.frame(
      Dataset = "Targeted",
      Category = "IGHD alleles",
      Expected = 10,
      Observed = sum(Segment == "IGHD" & Category %in% "Functional")),
    # "IGHJ: 56 IGHJ genes were found in the nine macaques including 9 unique
    # IGHJ genes and 11 alleles. Two alleles were identified for the IGHJ1 and
    # IGHJ3P genes. No IGHJ genes were identified that were not found in the
    # reference assembly. In no monkeys were two copies of IGHJ5 found, although
    # the reference does contain two copies, as noted above. With the exception
    # of M1, the single copy of IGHJ5 was functional. For M1, the entire IGHJ
    # cluster was assembled onto one contig, containing eight of the nine IGHJ
    # found in the reference assembly (4 functional) and one IGHJ5P pseudogene.
    # The complete IGHJ loci found in M0 and M1 were approximately 99.4%
    # identical (Figure S2 in Supplementary Material) apart from a 402 bp
    # indel."
    data.frame(
      Dataset = "Targeted",
      Category = "IGHJ alleles",
      Expected = 2,
      Observed = sum(Segment == "IGHJ" & Dataset == "Targeted")),
    # "IGHC: Owing to low coverage, constant region genes were assembled into
    # multiple contigs, typically distributing exons of the same gene into
    # different contigs. This circumstance makes it impossible to unambiguously
    # determine allelism. Nevertheless, the exons found on each of the
    # supplementary assemblies were identical with or very similar to exons on
    # each of the other supplementary assemblies and to exons on the reference
    # assembly (Table S5B in Supplementary Material). While the IGHD genes were
    # identical in all ten macaques, the other IGHC genes (IGHA, IGHE, IGHM,
    # IGHG) were highly diverse between macaques (Table 2)."
    # (TODO checks for IGHC)
    
    # "IGKV: 359 Functional/ORF IGKV genes were found in M1–M9, including 177
    # unique, functional genes, and 6 unique ORFs (Table S6A in Supplementary
    # Material, Figure 5). The IGKV allele library contains a total of 70
    # functional genes with 214 alleles. Three IGKV ORF’s with 6 alleles were
    # identified (Table 3)."
    data.frame(
      Dataset = "Both",
      Category = "Functional IGKV genes",
      Expected = 70,
      Observed = length(unique(subset(genbank_alleles, Segment == "IGKV" & Category %in% "Functional")$Gene))),
    data.frame(
      Dataset = "Both",
      Category = "Functional IGKV alleles",
      Expected = 214,
      Observed = nrow(subset(genbank_alleles, Segment == "IGKV" & Category %in% "Functional"))),
    data.frame(
      Dataset = "Both",
      Category = "ORF IGKV genes",
      Expected = 3,
      Observed = length(unique(subset(genbank_alleles, Segment == "IGKV" & Category %in% "ORF")$Gene))),
    data.frame(
      Dataset = "Both",
      Category = "ORF IGKV alleles",
      Expected = 6,
      Observed = nrow(subset(genbank_alleles, Segment == "IGKV" & Category %in% "ORF"))),
    # "IGKJ: The IGKJ cluster was assembled into one contig each for each of the
    # nine macaques. Each such cluster contained four functional genes and one
    # pseudogene (IGKJ5), with similarity among alleles uniformly above 99%. The
    # IGKJ clusters in two monkeys (M2, M8) are identical to that of the
    # reference."
    
    # "IGKC: The single IGKC was found in all animals in six unique allelic
    # forms. The assembly containing IGKC in four monkeys (M2, M4, M6, and M9)
    # are identical to that of the reference."
    data.frame(
      Dataset = "Both",
      Category = "IGKC genes",
      Expected = 1,
      Observed = length(unique(subset(genbank_alleles, Locus == "IGK" & Region == "constant")$Gene))),
    data.frame(
      Dataset = "Both",
      Category = "IGKC alleles",
      Expected = 6,
      Observed = nrow(subset(genbank_alleles, Locus == "IGK" & Region == "constant"))),
    # "IGLV: 291 genes were found, including 125 unique functional genes and 5
    # ORFs (Table S6B in Supplementary Material). The IGLV allele library
    # contains 50 functional genes with 137 alleles and 2 ORFs with 5 alleles
    # (Table 4, Figure 6). Similarity of genes within each family varied by
    # family."
    data.frame(
      Dataset = "Both",
      Category = "Functional IGLV genes",
      Expected = 50,
      Observed = length(unique(subset(genbank_alleles, Segment == "IGLV" & Category %in% "Functional")$Gene))),
    data.frame(
      Dataset = "Both",
      Category = "Functional IGLV alleles",
      Expected = 137,
      Observed = length(unique(subset(genbank_alleles, Segment == "IGLV" & Category %in% "Functional")$Allele)),
      Comments = "Figure 4 shows 140 functional alleles"),
    data.frame(
      Dataset = "Both",
      Category = "ORF IGLV genes",
      Expected = 2,
      Observed = length(unique(subset(genbank_alleles, Segment == "IGLV" & Category %in% "ORF")$Gene))),
    data.frame(
      Dataset = "Both",
      Category = "ORF IGLV alleles",
      Expected = 5,
      Observed = length(unique(subset(genbank_alleles, Segment == "IGLV" & Category %in% "ORF")$Allele))),
    # "IGLJ: 41 genes were found including 10 unique functional genes. IGL4 and
    # IGL5 had mutated recombination signal sequences and were NF. The IGLJ
    # allele library contains 7 IGLJ with 11 alleles."
    data.frame(
      Dataset = "Both",
      Category = "Total IGLJ genes",
      Expected = 7,
      Observed = length(unique(subset(genbank_alleles, Segment == "IGLJ")$Gene))),
    data.frame(
      Dataset = "Both",
      Category = "Total IGLJ alleles",
      Expected = 11,
      Observed = length(unique(subset(genbank_alleles, Segment == "IGLJ")$Allele))),
    # "IGLC: Very few IGLC genes were identified in M1–M9. 16 genes were found,
    # including 12 unique sequences. Of the 12 alleles identified, only 2
    # (IGLC2/IGLC3 and IGLC4P) of these were identical to the IGLC genes
    # identified in the reference assembly. The final IGLC allele library
    # consisted of 7 IGLJ genes (17 alleles)."
    data.frame(
      Dataset = "Both",
      Category = "Total IGLC genes",
      Expected = 7,
      Observed = length(unique(subset(genbank_alleles, Locus == "IGL" & Region == "constant")$Gene))),
    data.frame(
      Dataset = "Both",
      Category = "Total IGLC alleles",
      Expected = 17,
      Observed = length(unique(subset(genbank_alleles, Locus == "IGL" & Region == "constant")$Allele)))),
    function(df) {
      if (! "Comments" %in% colnames(df)) {
        df$Comments <- ""
      }
      df
    }
  )))
  checks <- rbind(wgs_checks, other_checks)
  checks
}


# Other -------------------------------------------------------------------


# A few GenBank entries have inconsistencies that I think I can resolve here.
# Should point them out at some point though.
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

# What was I using this for again?
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