
wildcard_constraints:
    acc="[A-Z]{2}[0-9]+"

with open("from-paper/accessions.txt") as f_in:
    ACCESSIONS = [line.strip() for line in f_in]

GBF = expand("from-genbank/{acc}.gbf", acc=ACCESSIONS)
FASTA = expand("from-genbank/{acc}.fasta", acc=ACCESSIONS)
GBFFASTA = expand("converted/{acc}.gbf.fasta", acc=ACCESSIONS)
GBFCSV = expand("converted/{acc}.gbf.csv", acc=ACCESSIONS)

SHEETS_URL = "https://docs.google.com/spreadsheets/d/e/2PACX-1vQC9Yh2l3ktopxK0idgGlSaOeo2Chnq15EoVro3wsmKHowVP1uVydyVJ_asCQe9Sfwot7_PcTNzaKGa/pub"
SHEETS_GIDS = {
    "fig1": "1857954478",
    "fig2": "380706608",
    "fig3": "1042301659",
    "fig4": "1508567143",
    "fig5": "782424017",
    "fig6": "1329373153",
    "suppsheet1": "974273021",
    "suppsheet2": "1621508193",
    "suppsheet3": "664446886",
    "table1": "655585230",
    "table3": "78785043",
    "table4": "1988184408"
    }
SHEETS = expand("from-paper/{sheet}.csv", sheet=SHEETS_GIDS.keys())

rule alleles_fasta:
    output: "output/alleles.fasta"
    input: "output/alleles.csv"
    shell: "python scripts/csv_to_fasta.py {input} {output} Allele SeqCDS"

rule merge_everything:
    output: "output/alleles.csv"
    input:
        from_gb="converted/all.csv",
        from_paper=SHEETS
    shell: "Rscript scripts/merge_everything.R"

rule convert_gbf_csv_combined:
    """Merge the per-accession CSV files with per-feature rows into one CSV."""
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

rule all_gbf_csv:
    input: GBFCSV

rule all_gbf_fasta:
    input: GBFFASTA

rule all_download_gbf:
    input: GBF

rule all_download_fasta:
    input: FASTA

rule all_sheets:
    input: SHEETS

rule convert_gbf_csv:
    """Convert a GBF file into a CSV with one row per feature."""
    output: "converted/{acc}.gbf.csv"
    input: "from-genbank/{acc}.gbf"
    shell: "python scripts/convert_gbf.py {input} {output} csv"

rule convert_gbf_fasta:
    """Convert a GBF file into a FASTA with one sequence per feature."""
    output: "converted/{acc}.gbf.fasta"
    input: "from-genbank/{acc}.gbf"
    shell: "python scripts/convert_gbf.py {input} /dev/stdout fasta | seqtk seq -l 0 > {output}"

rule download_gbf:
    """Download one GBF text file per GenBank accession.

    We can then convert the GBF into other formats like individual feature
    sequences in FASTA.
    """
    output: "from-genbank/{acc}.gbf"
    shell: "python scripts/download_genbank.py {wildcards.acc} gb > {output}"

rule download_fasta:
    """Download one FASTA per GenBank accession.

    Not actively using this since what I really want is to slice up each
    accession into one sequence per feature, to get at the individual
    genes/alleles.
    """
    output: "from-genbank/{acc}.fasta"
    shell: "python scripts/download_genbank.py {wildcards.acc} fasta > {output}"

rule download_sheet:
    output: "from-paper/{sheet}.csv"
    params:
        url=SHEETS_URL,
        gid=lambda w: SHEETS_GIDS[w.sheet]
    shell:
        """
            curl -L '{params.url}?gid={params.gid}&single=true&output=csv' > {output}
            dos2unix {output}
            echo >> {output}
        """
