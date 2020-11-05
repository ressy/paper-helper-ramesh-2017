
with open("accessions.txt") as f_in:
    ACCESSIONS = [line.strip() for line in f_in]

GBF = expand("from-genbank/{acc}.gbf", acc=ACCESSIONS)
FASTA = expand("from-genbank/{acc}.fasta", acc=ACCESSIONS)

rule all_download_gbf:
    input: GBF

rule all_download_fasta:
    input: FASTA

rule download_gbf:
    output: "from-genbank/{acc}.gbf"
    shell: "python download_genbank.py {wildcards.acc} gb > {output}"

rule download_fasta:
    output: "from-genbank/{acc}.fasta"
    shell: "python download_genbank.py {wildcards.acc} fasta > {output}"
