#!/usr/bin/env python3
import argparse
import gzip
import re
import subprocess
import sys

from functools import lru_cache

from Bio import Entrez
from ete3 import NCBITaxa


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


class LineageLookup:

    TAXLEVELS = {
        "kingdom": (0, "k"),
        "phylum": (1, "p"),
        "class": (2, "c"),
        "order": (3, "o"),
        "family": (4, "f"),
        "genus": (5, "g"),
        "species": (6, "s"),
    }

    def __init__(self):
        self.ncbi = NCBITaxa()

    @lru_cache(maxsize=1000)
    def get_lineage(self, taxid):


        def get_lineage_str(levels):
            return ";".join(
                f"{prefix}__{string}"
                for prefix, string in sorted(levels)
            )


        try:
            lineage = self.ncbi.get_lineage(taxid)
        except ValueError:
            return (taxid, "INVALID", "INVALID")

        ranks = self.ncbi.get_rank(lineage)
        rev_ranks = {v: k for k, v in ranks.items()}
        names = self.ncbi.get_taxid_translator(lineage)

        levels = []
        for tid, tname in zip(lineage, self.ncbi_translate_to_names(lineage)):
            rank = ranks.get(tid, "").replace("superkingdom", "kingdom")
            try:
                index, prefix = LineageLookup.TAXLEVELS.get(rank)
            except TypeError:
                continue
            levels.append((index, prefix, str(tid), tname))

        return get_lineage_str((index, prefix, tid) for index, prefix, tid, _ in levels), get_lineage_str((index, prefix, tname) for index, prefix, _, tname in levels)




#>>> ncbi.get_lineage(662479)
#[1, 131567, 2157, 28890, 2290931, 183963, 1644055, 1644056, 2251, 403181, 662479]
#>>> ncbi.get_rank(ncbi.get_lineage(662479))
#{1: 'no rank', 2157: 'superkingdom', 2251: 'genus', 28890: 'phylum', 131567: 'no rank', 183963: 'class', 403181: 'species', 662479: 'strain', 1644055: 'order', 1644056: 'family', 2290931: 'clade'}
#>>> ncbi.translate_to_names(ncbi.get_lineage(662479))
#['root', 'cellular organisms', 'Archaea', 'Euryarchaeota', 'Stenosarchaea group', 'Halobacteria', 'Haloferacales', 'Haloferacaceae', 'Haloferax', 'Haloferax mucosum', 'Haloferax mucosum ATCC BAA-1512']
#>>> ncbi.get_taxid_translator(ncbi.get_lineage(662479))
#{1: 'root', 2157: 'Archaea', 2251: 'Haloferax', 28890: 'Euryarchaeota', 131567: 'cellular organisms', 183963: 'Halobacteria', 403181: 'Haloferax mucosum', 662479: 'Haloferax mucosum ATCC BAA-1512', 1644055: 'Haloferacales', 1644056: 'Haloferacaceae', 2290931: 'Stenosarchaea group'}
#




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
    taxids = set()
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
            print(acc, mod_acc, acc_data, *lineage_lookup.get_lineage(acc_data["taxid"]))

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
