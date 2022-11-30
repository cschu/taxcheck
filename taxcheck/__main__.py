#!/usr/bin/env python3
import argparse
import gzip
import re
import subprocess
import sys

from taxcheck.lineage import LineageLookup



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
            yield aln[0], tuple(tags.get("YP").split(",")), tuple(tags.get("XA", "").split(";"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bamfile", type=str)
    
    args = ap.parse_args()

    cmd = ("samtools", "view", "-F", "0x4", args.bamfile)
    sam_proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    lineage_lookup = LineageLookup()

    for rname, yp_tag, xa_tag in extract_yp_reads_from_sam(sam_proc.stdout):
        print(rname, yp_tag, len(xa_tag))
        lineages = tuple(lineage_lookup.get_lineage(taxid) for taxid in yp_tag)
        for _, _, lineage in lineages:
            print(lineage)

    
    


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
