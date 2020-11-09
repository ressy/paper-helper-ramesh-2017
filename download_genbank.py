#!/usr/bin/env python
"""
Simple GenBank accession downloader.
"""

import sys
import time
from urllib.error import HTTPError
from Bio import Entrez

Entrez.email = "ancon@pennmedicine.upenn.edu"

def _get_gb_entry(**kwargs):
    args = {"db": "nucleotide"}
    args.update(kwargs)
    while True:
        try:
            handle = Entrez.efetch(**args)
        except HTTPError as error:
            sys.stderr.write(str(error))
            sys.stderr.write("\n")
            time.sleep(5)
        else:
            break
    return handle

def download_genbank(acc, rettype, out_handle=sys.stdout):
    """Download a file for a single genbank accession.

    acc: GenBank accession
    rettype: "fasta" or "gb"
    out_handle: an open file handle to write output to
    """
    handle = _get_gb_entry(id=acc, rettype=rettype, retmode="text")
    out_handle.write(handle.read().strip() + "\n")

if __name__ == "__main__":
    download_genbank(sys.argv[1], sys.argv[2])
