# Ramesh paper's antibody alleles versus SONAR's copy

Where is `IGHV2-ABO*01` in SONAR's germDB files?  Ryan found that that was a
close match for RM6561's DH1021 lineage, but it's not in the IgDiscover results
and in turn is not in the initial Ramesh database supplied by SONAR even though
it's visible in GenBank. I should compare the original Ramesh files in GenBank
with what's in SONAR's copy, in a systematic way.

<https://www.frontiersin.org/articles/10.3389/fimmu.2017.01407/full>

accessions MF989451 to MF989952

First glance: downloading a FASTA per accession does *not* reveal
`IGHV2-ABO*01` despite that allele being shown in the GenBank entry for
MF989491.1.  It looks like that's one of the combo entries labeled "locus" and
there isn't a separate "gene" one for it.  Is that why it's not in Chaim's set,
maybe?
