#!/usr/bin/env python
"""
Convert GBF features to CSV.

Note that this is not the same as just using the SeqRecord provided by SeqIO's
GenBank parser.  That produces one big sequence for the whole record.  Instead
this writes each feature as a separate CSV row.
"""

import sys
import csv
import itertools
from Bio import SeqIO
from Bio.SeqFeature import BeforePosition, AfterPosition, ExactPosition

def get_gbf_attrs(gbf):
    """Extract attributes of an entire GBF record."""
    attrs = {}
    attrs.update(gbf.annotations)
    del attrs["structured_comment"]
    del attrs["references"]
    attrs["seqlen"] = len(gbf.seq)
    attrs["dbxrefs"] = gbf.dbxrefs
    attrs["description"] = gbf.description
    attrs["id"] = gbf.id
    attrs["name"] = gbf.name
    attrs["accessions"] = ";".join(attrs["accessions"])
    attrs["taxonomy"] = ";".join(attrs["taxonomy"])
    attrs["keywords"] = ";".join([k for k in attrs["keywords"] if k])
    attrs = {key: val for key, val in attrs.items() if val}
    return attrs

def get_feature_attrs(feature, gbf):
    """Extract attributes for one GBF feature."""
    # lookup table mapping the feature start/end position class to a
    # one-character label.  (We'll get a KeyError if the location objects
    # aren't instances of one of these handled cases but these should account
    # for everything we'll see here, I think.)
    location_lut = {
        ExactPosition: "=",
        BeforePosition: "<",
        AfterPosition: ">"}
    feature_pairs = {
        "id": feature.id,
        "ref": feature.ref,
        "ref_db": feature.ref_db,
        "strand": feature.strand,
        "type": feature.type,
        "pos_start": feature.location.start.position + 1,
        "pos_end": feature.location.end.position,
        "pos_start_type": location_lut[type(feature.location.start)],
        "pos_end_type": location_lut[type(feature.location.end)],
        "pos_strand": feature.location.strand,
        "seq": feature.extract(gbf.seq)}
    for key in feature.qualifiers:
        feature_pairs["qualifier_" + key] = ";".join(feature.qualifiers[key])
    if feature_pairs["id"] == "<unknown id>":
        feature_pairs["id"] = None
    feature_pairs = {key: val for key, val in feature_pairs.items() if val is not None}
    return feature_pairs

def read_gbf_as_table(fp_in):
    """Convert a GBF into a list of dictionaries, one per feature."""
    rows = []
    with open(fp_in) as f_in:
        for gbf in SeqIO.parse(f_in, "gb"):
            gbf_pairs = get_gbf_attrs(gbf)
            gbf_pairs = {"gbf_" + key: val for key, val in gbf_pairs.items()}
            feature_src_attrs = {}
            for feature in gbf.features:
                feature_pairs = get_feature_attrs(feature, gbf)
                feature_pairs = {"feature_" + key: val for key, val in feature_pairs.items()}
                # the source feature is special.  stash those attributes and
                # use as defaults for the other features.
                if feature.type == "source":
                    del feature_pairs["feature_seq"]
                    del feature_pairs["feature_type"]
                    feature_src_attrs = {
                        k.replace("feature", "feature_source"): v for k, v in feature_pairs.items()}
                else:
                    row = feature_src_attrs.copy()
                    row.update(feature_pairs)
                    row.update(gbf_pairs)
                    rows.append(row)
    return rows

def convert_gbf(fp_in, fp_out):
    """Convert a GBF file to CSV, with one row per feature."""
    rows = read_gbf_as_table(fp_in)
    with open(fp_out, "wt") as f_out:
        fieldnames = [row.keys() for row in rows]
        fieldnames = sorted(list(set(itertools.chain(*fieldnames))))
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            for key in fieldnames:
                if key not in row:
                    row[key] = ""
            writer.writerow(row)

if __name__ == "__main__":
    convert_gbf(sys.argv[1], sys.argv[2])
