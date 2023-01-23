#!/usr/bin/env python3
import argparse
import gzip
import json
import os
import re
import subprocess
import sys

from collections import Counter

from Bio import Entrez

from taxcheck.lineage import LineageLookup, LineageFactory, Lineage
from taxcheck.ncbi import ncbi_tax_lookup


def get_lines_from_chunks(_in, bufsize=400000000):
    tail = ""
    while True:
        chunk = "".join((tail, _in.read(bufsize).decode()))
        if not chunk:
            break
        chunk = chunk.split("\n")
        *chunk, tail = chunk
        for line in chunk:
            if line:
                yield line
    if tail:
        yield tail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fastafile", type=str)
    ap.add_argument("email", type=str)
    ap.add_argument("--debug", action="store_true")
    ap.add_argument("--ncbi-chunksize", type=int, default=400)
    args = ap.parse_args()

    if args.ncbi_chunksize < 10:
        raise ValueError("NCBI chunk size needs to be at least 10.")

    print("Extracting NCBI ids from fasta...", file=sys.stderr)
    with open(args.fastafile, "rb") as _in:
        refs = {
            line.strip()[1:].split(" ")[0]
            for line in get_lines_from_chunks(_in)
            if line[0] == ">"
        }
    
        print(f"Looking up {len(refs)} taxonomy ids...", file=sys.stderr) 
        ncbi_lookup = ncbi_tax_lookup(args.email, list(refs), chunksize=args.ncbi_chunksize)

        # for v in ncbi_lookup.values():
        #    print(v["accession"], v["id"], v["taxid"], sep="\t")


if __name__ == "__main__":
    main()
