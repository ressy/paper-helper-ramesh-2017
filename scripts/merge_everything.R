#!/usr/bin/env Rscript
source("scripts/analysis.R")
genbank_entries <- read.csv("converted/all.csv", stringsAsFactors = FALSE)
genbank_alleles <- parse_genbank_entries(genbank_entries)
paper <- load_paper("from-paper")
metadata <- parse_paper(paper)
genbank_alleles <- merge_metadata(genbank_alleles, metadata)
genbank_alleles <- collapse_fields(genbank_alleles)
genbank_alleles <- finalize_table(genbank_alleles)
checks <- make_check_tables(genbank_alleles)
for (tblname in names(checks)) {
  write.csv(
    checks[[tblname]],
    file.path("checks", paste0(tblname, ".csv")),
    row.names = FALSE,
    quote = FALSE,
    na = "")
}
write.csv(genbank_alleles, "output/alleles.csv", row.names=FALSE, na = "")
