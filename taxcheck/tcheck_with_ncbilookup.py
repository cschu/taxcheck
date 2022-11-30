#!/usr/bin/env python3
import argparse
import gzip
import re
import subprocess
import sys

from Bio import Entrez

from lineage import LineageLookup



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


def extract_read_ref_from_sam(stream):
    read2ref = {}
    for aln in get_lines_from_chunks(stream):
        aln = aln.strip().split("\t")
        flag = int(aln[1])
        if flag & 0x1 and not (aln[2] == aln[6] or flag & 0x8):
            continue
        tags = dict(item.split(":")[0::2] for item in aln[11:])
        if not tags.get("XA"):
            read2ref.setdefault(aln[2], set()).add(aln[0])
    return read2ref






def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bamfile", type=str)
    ap.add_argument("email", type=str)
    args = ap.parse_args()

    Entrez.email = args.email

    cmd = ("samtools", "view", "-F", "0xf04", args.bamfile)
    sam_proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    read2ref = extract_read_ref_from_sam(sam_proc.stdout)


    ncbi_lookup = {}
    lineage_lookup = LineageLookup()

    accessions = tuple(read2ref)
    chunksize = 20
    for start in range(0, len(read2ref), chunksize):
        ids = accessions[start:start + chunksize]
        try:
            epost_results = Entrez.read(Entrez.epost(db="nucleotide", id=",".join(ids), format="acc"))
        except RuntimeError as err:
            print("Problem with IDs:", *ids, sep="\n", file=sys.stderr)
            continue
        efetch_handle = Entrez.efetch(db="nucleotide", rettype="docsum", retmode="xml", start=start, retmax=chunksize, webenv=epost_results["WebEnv"], query_key=epost_results["QueryKey"], idtype="acc")
        data = Entrez.read(efetch_handle)

        ncbi_lookup = {
            item["AccessionVersion"]: {
                "accession": item["AccessionVersion"],
                "id": item["Id"],
                "taxid": int(item["TaxId"]),
            }
            for item in data
        }

        for acc in ids:
            mod_acc = re.sub(r"^ref\|", "", acc).strip("|")
            acc_data = ncbi_lookup.get(mod_acc)
            taxname, name_lineage, taxid_lineage = lineage_lookup.get_lineage(acc_data["taxid"])
            # print(acc, mod_acc, acc_data, name_lineage, taxid_lineage)
            
            for read in read2ref[acc]:
                print(read, acc, mod_acc, acc_data["taxid"], taxname, name_lineage, taxid_lineage, sep="\t")
            

#[
#    {
#        'Item': [], 'Id': '1160384904', 'Caption': 'NZ_MZGV01000142',
#        'Title': 'Clostridium oryzae strain DSM 28571 CLORY_contig000142, whole genome shotgun sequence',
#        'Extra': 'gi|1160384904|ref|NZ_MZGV01000142.1||gnl|WGS:NZ_MZGV01|CLORY_contig000142[1160384904]',
#        'Gi': IntegerElement(1160384904, attributes={}), 'CreateDate': '2017/03/17', 'UpdateDate': '2022/04/01',
#        'Flags': IntegerElement(544, attributes={}), 'TaxId': IntegerElement(1450648, attributes={}),
#        'Length': IntegerElement(2352, attributes={}), 'Status': 'live', 'ReplacedBy': '', 'Comment': '  ', 'AccessionVersion': 'NZ_MZGV01000142.1'
#    },
#    {
#        'Item': [], 'Id': '1180753414', 'Caption': 'NZ_FWXW01000002',
#        'Title': 'Papillibacter cinnamivorans DSM 12816, whole genome shotgun sequence',
#        'Extra': 'gi|1180753414|ref|NZ_FWXW01000002.1|[1180753414]', 'Gi': IntegerElement(1180753414, attributes={}), 'CreateDate': '2017/04/08', 'UpdateDate': '2022/04/01',
#        'Flags': IntegerElement(544, attributes={}), 'TaxId': IntegerElement(1122930, attributes={}),
#        'Length': IntegerElement(423394, attributes={}), 'Status': 'live', 'ReplacedBy': '', 'Comment': '  ', 'AccessionVersion': 'NZ_FWXW01000002.1'
#    }
#]


if __name__ == "__main__":
    main()
