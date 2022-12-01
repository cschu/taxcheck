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
    
    args = ap.parse_args()

    cmd = ("samtools", "view", "-F", "0x4", args.bamfile)
    sam_proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    lineage_lookup = LineageLookup()
    lfactory = LineageFactory()

    read2ref = {}
    refs = set()
    
    print("Reading bam file...")
    for rname, ref, mapq, yp_tag, xa_tag in extract_yp_reads_from_sam(sam_proc.stdout):
        read2ref[rname] = {
            "ref": ref,
            "mapq": mapq,
            "yp": yp_tag,
            "xa": xa_tag,
        }
        refs.update((ref, ) + (tuple(item[0] for item in xa_tag) if xa_tag else tuple()))
        # print(rname, yp_tag, len(xa_tag))
        # lineages = tuple(lineage_lookup.get_lineage(taxid) for taxid in yp_tag)
        # for _, _, lineage in lineages:
        #    print(lineage)
        if len(read2ref) > 10:
            break

    with open("refs.txt", "wt") as _out:
        print(*sorted(refs), sep="\n", file=_out)
       
    print(f"Looking up {len(refs)} taxonomy ids...") 
    cached_ncbi = "ncbi_lut.json"
    if os.path.isfile(cached_ncbi):
        with open(cached_ncbi) as _in:
            ncbi_lookup = json.load(_in)
    else:
        ncbi_lookup = ncbi_tax_lookup(args.email, list(refs), chunksize=400)
        with open(cached_ncbi, "wt") as _out:
            json.dump(ncbi_lookup, _out)

    print("Annotating reads...")
    for rname, aln_data in read2ref.items():

        # targets = (aln_data["ref"],) + tuple(ref for ref, _ in aln_data["xa"]) if aln_data["xa"] else tuple()
        # print(targets)
        # lcount = Counter
        print(rname, aln_data)
        rrefs = ((aln_data["ref"],) + (tuple(r for r, _ in aln_data["xa"]) if aln_data["xa"] else tuple()))
        print("RREFS", rrefs)
        print("CHECK", aln_data["ref"] in ncbi_lookup)

        lcount = Counter()
        for ref in rrefs:
            print(ref, ncbi_lookup.get(ref))
            lcount[ncbi_lookup.get(ref, {}).get("taxid", -1)] += 1
        print("LCOUNT", lcount)
        print("LCOUNT HAS -1", -1 in lcount)
        #break 

        #lcount = Counter(
        #    ncbi_lookup.get(ref, ncbi_tax_lookup(args.email, [ref])).get("taxid", -1)
        #    for ref in ((aln_data["ref"],) + (tuple(r for r, _ in aln_data["xa"]) if aln_data["xa"] else tuple()))
        #)
        #break
        lineage_data = {
            taxid: {
                "lineage": lfactory.generate_lineage(taxid),
                "count": lcount[taxid], 
                #sum(
                #    1 
                #    for ref, _ in aln_data["xa"] + (aln_data["ref"], aln_data["mapq"])
                #    if ncbi_lookup.get(ref, {}).get("taxid", -1) == taxid
                #),
            }
            for taxid in aln_data["yp"]
        }
        print(lineage_data)
        lineages = tuple(tuple(v.values()) for v in lineage_data.values())  #  for item in lineage_data.values()] 
        print(lineages)
        print(*(f"{l.get_string()}: {c}" for l, c in lineages), sep="\n")
        if len(lineages) == 1:
            consensus_lineage = lineages[0][0]
        else:
            for level in range(5, -1, -1):
                c = Counter()
                for lineage, count in lineages:
                    c.update((lineage.levels[level]["id"],) * count)
                print(c)
                if c:
                    top_taxid, top_count = c.most_common()[0]
                    if top_count / len(c) > 0.75:
                        consensus_lineage = lfactory.generate_lineage(top_taxid)
                        break
        print(rname, consensus_lineage.get_string(), consensus_lineage.get_string(show_names=False), sep="\t")



    


    # accessions = tuple(read2ref)
    # chunksize = 20
    # for start in range(0, len(read2ref), chunksize):
    #     ids = accessions[start:start + chunksize]
    #     try:
    #         epost_results = Entrez.read(Entrez.epost(db="nucleotide", id=",".join(ids), format="acc"))
    #     except RuntimeError as err:
    #         print("Problem with IDs:", *ids, sep="\n", file=sys.stderr)
    #         continue
    #     efetch_handle = Entrez.efetch(db="nucleotide", rettype="docsum", retmode="xml", start=start, retmax=chunksize, webenv=epost_results["WebEnv"], query_key=epost_results["QueryKey"], idtype="acc")
    #     data = Entrez.read(efetch_handle)

    #     ncbi_lookup = {
    #         item["AccessionVersion"]: {
    #             "accession": item["AccessionVersion"],
    #             "id": item["Id"],
    #             "taxid": int(item["TaxId"]),
    #         }
    #         for item in data
    #     }

    #     for acc in ids:
    #         mod_acc = re.sub(r"^ref\|", "", acc).strip("|")
    #         acc_data = ncbi_lookup.get(mod_acc)
    #         taxname, name_lineage, taxid_lineage = lineage_lookup.get_lineage(acc_data["taxid"])
    #         # print(acc, mod_acc, acc_data, name_lineage, taxid_lineage)
            
    #         for read in read2ref[acc]:
    #             print(read, acc, mod_acc, acc_data["taxid"], taxname, name_lineage, taxid_lineage, sep="\t")
            


if __name__ == "__main__":
    main()
