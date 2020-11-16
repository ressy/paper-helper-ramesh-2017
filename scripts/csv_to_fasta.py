#!/usr/bin/env python
"""Extract sequences as FASTA from two columns of CSV."""

import sys
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def csv_to_fasta(fp_csv_in, fp_fasta_out, seqid_col, seq_col):
    """Extract sequences as FASTA from two columns of CSV."""
    with open(fp_csv_in) as f_in, open(fp_fasta_out, "wt") as f_out:
        reader = csv.DictReader(f_in)
        for row in reader:
            record = SeqRecord(Seq(row[seq_col]), id=row[seqid_col], description="")
            SeqIO.write(record, f_out, "fasta")

if __name__ == "__main__":
    csv_to_fasta(*sys.argv[1:])
