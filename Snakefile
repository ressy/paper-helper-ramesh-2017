
wildcard_constraints:
    acc="[A-Z]{2}[0-9]+"

with open("accessions.txt") as f_in:
    ACCESSIONS = [line.strip() for line in f_in]

GBF = expand("from-genbank/{acc}.gbf", acc=ACCESSIONS)
FASTA = expand("from-genbank/{acc}.fasta", acc=ACCESSIONS)
GBFFASTA = expand("converted/{acc}.gbf.fasta", acc=ACCESSIONS)
GBFCSV = expand("converted/{acc}.gbf.csv", acc=ACCESSIONS)

rule all_gbf_csv:
    input: GBFCSV

rule all_gbf_fasta:
    input: GBFFASTA

rule all_download_gbf:
    input: GBF

rule all_download_fasta:
    input: FASTA

rule convert_gbf_csv_combined:
    output: "converted/all.csv"
    input: GBFCSV
    run:
        import csv
        rows = []
        csv.field_size_limit(10000000)
        for fp in input:
            with open(fp) as f_in:
                reader = csv.DictReader(f_in)
                rows.extend(list(reader))
        with open(output[0], "wt") as f_out:
            writer = csv.DictWriter(f_out, fieldnames=rows[0].keys(), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

rule convert_gbf_csv:
    output: "converted/{acc}.gbf.csv"
    input: "from-genbank/{acc}.gbf"
    shell: "python convert_gbf.py {input} {output} csv"

rule convert_gbf_fasta:
    output: "converted/{acc}.gbf.fasta"
    input: "from-genbank/{acc}.gbf"
    shell: "python convert_gbf.py {input} /dev/stdout fasta | seqtk seq -l 0 > {output}"

rule download_gbf:
    output: "from-genbank/{acc}.gbf"
    shell: "python download_genbank.py {wildcards.acc} gb > {output}"

rule download_fasta:
    output: "from-genbank/{acc}.fasta"
    shell: "python download_genbank.py {wildcards.acc} fasta > {output}"
