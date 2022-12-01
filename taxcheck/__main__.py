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
            yield line
    if tail:
        yield tail


def extract_yp_reads_from_sam(stream):    
    for aln in get_lines_from_chunks(stream):
        aln = aln.strip().split("\t")
        tags = dict(item.split(":")[0::2] for item in aln[11:])
        if tags.get("YP"):
            xa_tag = tags.get("XA")
            if xa_tag is None:
                xa_tag = tuple()
            else:
                xa_tag = tuple(item.split(",")[0::3] for item in xa_tag.strip(";").split(";"))
            yield aln[0], aln[2], aln[4], tuple(map(int, tags.get("YP").split(","))), xa_tag


def extract_refs_from_xa_tag(xa_tag):
    # XA:Z:NZ_BAHC01000194.1,-269,60S30M,0;NZ_BAFX01000157.1,-297,60S29M1S,0;NZ_MOSG01000117.1,-656,59S25M6S,0;NZ_AOON01000169.1,+2,25M65S,0;
    return (item.split(",")[0::3] for item in xa_tag[5:].split(";"))    


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bamfile", type=str)
    ap.add_argument("email", type=str)
    ap.add_argument("--debug", action="store_true")
    ap.add_argument("--lineage-cutoff", type=float, default=0.75)
    ap.add_argument("--species-cutoff", type=float, default=0.99)
    args = ap.parse_args()

    cmd = ("samtools", "view", "-F", "0x4", args.bamfile)
    sam_proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    lfactory = LineageFactory()
    read2ref = {}
    refs = set()
    
    print("Reading bam file...", file=sys.stderr)
    for rname, ref, mapq, yp_tag, xa_tag in extract_yp_reads_from_sam(sam_proc.stdout):
        read2ref[rname] = {
            "ref": ref,
            "mapq": mapq,
            "yp": yp_tag,
            "xa": xa_tag,
        }
        refs.update((ref, ) + (tuple(item[0] for item in xa_tag) if xa_tag else tuple()))
        if args.debug and len(read2ref) > 10:
            break

    print(f"Looking up {len(refs)} taxonomy ids...", file=sys.stderr) 
    cached_ncbi = "ncbi_lut.json"
    if args.debug and os.path.isfile(cached_ncbi):
        with open(cached_ncbi) as _in:
            ncbi_lookup = json.load(_in)
    else:
        ncbi_lookup = ncbi_tax_lookup(args.email, list(refs), chunksize=400)
        with open(cached_ncbi, "wt") as _out:
            json.dump(ncbi_lookup, _out)

    print("Annotating reads...", file=sys.stderr)
    for rname, aln_data in read2ref.items():

        rrefs = ((aln_data["ref"],) + (tuple(r for r, _ in aln_data["xa"]) if aln_data["xa"] else tuple()))
        
        if args.debug:
            print(rname, aln_data, file=sys.stderr)
            print("RREFS", rrefs, file=sys.stderr)
            print("CHECK", aln_data["ref"] in ncbi_lookup, file=sys.stderr)

        # count how often a read aligns against a refseq from a specific taxonomy
        lcount = Counter()
        for ref in rrefs:
            if args.debug:
                print(ref, ncbi_lookup.get(ref), file=sys.stderr)
            lcount[ncbi_lookup.get(ref, {}).get("taxid", -1)] += 1
        
        if args.debug:
            print("LCOUNT", lcount, file=sys.stderr)
            print("LCOUNT HAS -1", -1 in lcount, file=sys.stderr)

        # generate lineages from each reference and store the alignment counts along with it  
        lineages2 = {
            taxid: {
                "lineage": lfactory.generate_lineage(taxid),
                "count": lcount[taxid], 
            }
            for taxid in aln_data["yp"]
        }

        if args.debug:
            print(*lineages2, sep="\n", file=sys.stderr)            
            print(*(f"{l.get_string()}: {c}" for l, c in lineages2.values()), sep="\n", file=sys.stderr)

        # determine a consensus line based on fraction of alignments
        if len(lineages2) == 1:
            consensus_lineage = tuple(lineages2.values())[0]["lineage"]
        else:
            # iterate from species -> kingdom
            for level in range(Lineage.TAXLEVELS["species"][0], -1, -1):
                tax_counter = Counter()
                # for each level count the taxids (multiplied by number of alignments)
                for lineage in lineages2.values():
                    tax_counter.update((lineage["lineage"].levels[level]["id"],) * lineage["count"])
                if args.debug:
                    print("LEVEL", level, tax_counter, file=sys.stderr)
                # then check if there's a consensus (based on fractional representation by the alignments)
                if tax_counter:
                    top_taxid, top_count = tax_counter.most_common()[0]
                    cutoff = args.species_cutoff if level == Lineage.TAXLEVELS["species"][0] else args.lineage_cutoff
                    if top_count / sum(tax_counter.values()) > cutoff:
                        consensus_lineage = lfactory.generate_lineage(top_taxid)
                        break    
                    
        print(rname, consensus_lineage.get_string(), consensus_lineage.get_string(show_names=False), sep="\t")


            


if __name__ == "__main__":
    main()
