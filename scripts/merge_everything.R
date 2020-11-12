#!/usr/bin/env Rscript
source("scripts/analysis.R")
genbank <- read.csv("converted/all.csv", stringsAsFactors = FALSE)
genbank_alleles <- parse_genbank_genes(genbank)
paper <- load_paper()
metadata <- list(
  genes = parse_genes(paper),
  scaffolds = parse_scaffolds(paper)
)
genbank_alleles <- merge_metadata(genbank_alleles, metadata)
write.csv(genbank_alleles, "output/alleles.csv", row.names=FALSE)
