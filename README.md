# Data Gathering from Ramesh 2017

Scripts to merge GenBank-deposited antibody sequences and annotations with
paper-supplied information (see below) into a single CSV table and FASTA, for
as many alleles as possible. See [output/alleles.csv](output/alleles.csv) and
[output/alleles.fasta](output/alleles.fasta) for final output (all else here is
an implementation detail).  The CSV file contains separate columns for genomic,
CDS, and AA sequences where available, and the FASTA file contains CDS
sequences.  **NB: The CSV file is still misssing some information, and shows
small inconsistencies compared to the paper I couldn't resolve here.**  I'll
make updates here, if I find any, to resolve these issues.

Ramesh A, Darko S, Hua A, Overman G, Ransier A, Francica JR, Trama A, Tomaras GD, Haynes BF, Douek DC and Kepler TB (2017) Structure and Diversity of the Rhesus Macaque Immunoglobulin Loci through Multiple De Novo Genome Assemblies. Front. Immunol. 8:1407. doi: [10.3389/fimmu.2017.01407](https://doi.org/10.3389/fimmu.2017.01407)
