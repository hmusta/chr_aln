#!/usr/bin/env python3

_COMPLEMENT = str.maketrans(
    "ACGTNacgtnRYSWKMBDHVryswkmbdhv",
    "TGCANtgcanYRSWMKVHDByrswmkvhdb",
)

def read_fasta(fname):
    with open(fname, 'r') as f:
        header = f.readline().rstrip().split()[0][1:]
        seq = "".join(line.rstrip().upper() for line in f.readlines())
    return header,seq
