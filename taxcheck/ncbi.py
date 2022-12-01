import sys

from Bio import Entrez


def ncbi_tax_lookup(email, accessions, chunksize=20):
    Entrez.email = email
    ncbi_lookup = {}

    for start in range(0, len(accessions), chunksize):
        ids = accessions[start:start + chunksize]
        if len(accessions) == 1:
            print(f"Looking up {ids[0]}...", end="")
        else:
            print(f"Looking up id block {start}-{min(len(accessions) - 1, start + chunksize - 1)}...", end="")
        try:
            epost_results = Entrez.read(Entrez.epost(db="nucleotide", id=",".join(ids), format="acc"))
        except RuntimeError as err:
            print("Problem with IDs:", *ids, sep="\n", file=sys.stderr)
            continue
        efetch_handle = Entrez.efetch(
            db="nucleotide", rettype="docsum", retmode="xml",
            idtype="acc", start=start, retmax=chunksize,
            webenv=epost_results["WebEnv"], query_key=epost_results["QueryKey"]
        )
        data = Entrez.read(efetch_handle)
        print(f"{len(data)} entries received.")

        for item in data:
            ncbi_lookup[item["AccessionVersion"]] = {
                "accession": item["AccessionVersion"],
                "id": item["Id"],
                "taxid": int(item["TaxId"]),
            }
    return ncbi_lookup
