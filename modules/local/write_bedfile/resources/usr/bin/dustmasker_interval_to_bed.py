#! /usr/bin/env python3

from csv import reader, writer
import sys


def get_id(s):
    # adapted from Biopython SeqIO fasta parser
    return s[1:].split(None, 1)[0]


r = reader(sys.stdin, delimiter="\t")
w = writer(sys.stdout, delimiter="\t")
for row in r:
    row[0] = get_id(
        row[0]
    )  # only keep the Accession number (trim everything after first space)
    row[2] = int(row[2]) + 1  # bed file 3rd field is 1-based
    w.writerow(row)