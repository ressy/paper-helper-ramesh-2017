#!/usr/bin/env Rscript
source("scripts/analysis.R")
genbank_entries <- read.csv("converted/all.csv", stringsAsFactors = FALSE)
genbank_alleles <- parse_genbank_entries(genbank_entries)
metadata <- parse_paper("from-paper")
genbank_alleles <- merge_metadata(genbank_alleles, metadata)
genbank_alleles <- collapse_fields(genbank_alleles)
genbank_alleles <- finalize_table(genbank_alleles)
write.csv(genbank_alleles, "output/alleles.csv", row.names=FALSE, na = "")