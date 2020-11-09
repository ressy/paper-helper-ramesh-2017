# Ramesh paper's antibody alleles versus SONAR's copy

Where is `IGHV2-ABO*01` in SONAR's germDB files?  Ryan found that that was a
close match for RM6561's DH1021 lineage, but it's not in the IgDiscover results
and in turn is not in the initial Ramesh database supplied by SONAR even though
it's visible in GenBank.  The Ramesh paper:

<https://www.frontiersin.org/articles/10.3389/fimmu.2017.01407/full>

The Snakemake rules download GBF files for MF989451 to MF989952 and extract
individual FASTA sequences for each feature, and the `analysis.R` file parses
out individual columns from the allele names and description text.  When
filtering to just complete sequences (Partial == FALSE) the numbers match
what's in the paper, but still not what's in SONAR's copy.  It doesn't look
like a straightforward pattern to me, like sequences missing from the
multi-gene accessions for example.  There are also some sequences in SONAR's
copy that don't occur in the GenBank copy, specifically ones with an extra
suffix on the name (like `IGKV3-ADU-S*01`) or possibly also with an "ORF"
prefix like `ORF_IGKV1-AES-X*01`.  Some of that might be explained by removing
those extra bits from the names but I'm not sure what they signify. 
