#!/usr/bin/env python3
import argparse
import gzip
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
            xa_tag = tuple(item.split(",")[0::3] for item in tags.get("XA", "").strip(";").split(";"))
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
    
    for rname, ref, mapq, yp_tag, xa_tag in extract_yp_reads_from_sam(sam_proc.stdout):
        read2ref[rname] = {
            "ref": ref,
            "mapq": mapq,
            "yp": yp_tag,
            "xa": xa_tag,
        }
        refs.update((ref, ) + tuple(item[0] for item in xa_tag))
        # print(rname, yp_tag, len(xa_tag))
        # lineages = tuple(lineage_lookup.get_lineage(taxid) for taxid in yp_tag)
        # for _, _, lineage in lineages:
        #    print(lineage)
        
    ncbi_lookup = ncbi_tax_lookup(args.email, list(refs))

    for rname, aln_data in read2ref.items():
        lineage_data = {
            taxid: {
                "lineage": lfactory.generate_lineage(taxid),
                "count": sum(
                    1 
                    for ref, _ in aln_data["xa"]
                    if ncbi_lookup.get(ref, {}).get("taxid", -1) == taxid
                ),
            }
            for taxid in aln_data["yp"]
        }
        lineages = [item["lineage"] for item in lineage_data.values()] 
        if len(lineages) == 1:
            consensus_lineage = lineages[0]
        else:
            for level in range(6, -1, -1):
                c = Counter(
                    lineage.levels[level]["id"] for lineage in lineages
                )
                if c:
                    top_taxid, top_count = c.most_common()[0]
                    if top_count / len(c) > 0.75:
                        consensus_lineage = lfactory.generate_lineage(top_taxid)
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
