import sys

from Bio import Entrez



def lookup_accessions(email, accessions, chunksize=20, db="nucleotide", trials=2):
    ncbi_lookup = {}

    try:
        epost_results = Entrez.read(Entrez.epost(db=db, id=",".join(accessions), format="acc"))
    except RuntimeError as err:
        if len(accessions) == 1:
            print(f"Cannot obtain lookup for id={accessions[0]}.", file=sys.stderr)
            return {accessions[0]: {"accession": accessions[0], "id": None, "taxid": None}}
        # print("Splitting accession block in half...", file=sys.stderr)
        accessions = sorted(accessions)
        n = len(accessions)
        block1, block2 = accessions[:n//2], accessions[n//2:]
        lookup1, lookup2 = lookup_accessions(email, block1, chunksize=chunksize, db=db), None
        if lookup1 is not None:
            ncbi_lookup.update(lookup1)
        if block2:
            lookup2 = lookup_accessions(email, block2, chunksize=chunksize, db=db)
            if lookup2 is not None:
                ncbi_lookup.update(lookup2)           
    else:
        efetch_handle = Entrez.efetch(
            db=db, rettype="docsum", retmode="xml", idtype="acc", retmax=chunksize,
            webenv=epost_results["WebEnv"], query_key=epost_results["QueryKey"],
        )
        try:
            data = Entrez.read(efetch_handle)
        except RuntimeError as err:
            
            if trials == 0:
                print(f"Caught runtime error when reading efetch:\n{err}\Lost accessions:", file=sys.stderr)
                print(",".join(accessions), file=sys.stderr)

                ncbi_lookup.update({
                    acc: {"accession": acc, "id": None, "taxid": None}
                    for acc in accessions
                })

            else:
                lookup = lookup_accessions(email, accessions, chunksize=chunksize, db=db, trials=trials - 1)
                if lookup is not None:
                    ncbi_lookup.update(lookup)
            
        else:
            print(f"Received {len(data)} entries.", file=sys.stderr)

            ncbi_lookup.update({
                item["AccessionVersion"]: {
                    "accession": item["AccessionVersion"],
                    "id": item["Id"],
                    "taxid": int(item["TaxId"]),
                }
                for item in data
            })

    return ncbi_lookup
        
             



def ncbi_tax_lookup(email, accessions, chunksize=20, db="nucleotide"):
    Entrez.email = email
    ncbi_lookup = {}

    for start in range(0, len(accessions), chunksize):
        ids = accessions[start:start + chunksize]
        if len(accessions) == 1:
            print(f"Looking up single id {ids[0]}...", file=sys.stderr)
        else:
            real_chunksize = min(len(accessions) - 1, start + chunksize - 1)
            print(f"Looking up id block {start}-{real_chunksize} of {len(accessions)} ({(real_chunksize) / len(accessions) * 100:.3f}%)", file=sys.stderr)

        lookup = lookup_accessions(email, accessions[start:start + chunksize], chunksize=chunksize, db=db)
        if lookup is not None:
            ncbi_lookup.update(lookup)
            for v in lookup.values():
                print(v["accession"], v["id"], v["taxid"], sep="\t")
            sys.stdout.flush()

        """
        try:
            epost_results = Entrez.read(Entrez.epost(db="nucleotide", id=",".join(ids), format="acc"))
        except RuntimeError as err:
            print("Problem with IDs:", "\n", file=sys.stderr)
            print(err, file=sys.stderr)
            print(*ids, sep=",", file=sys.stderr)
            raise RuntimeError from err
        efetch_handle = Entrez.efetch(
            db="nucleotide", rettype="docsum", retmode="xml",
            idtype="acc", start=start, retmax=chunksize,
            webenv=epost_results["WebEnv"], query_key=epost_results["QueryKey"]
        )
        data = Entrez.read(efetch_handle)
        print(f"{len(data)} entries received.", file=sys.stderr)

        for item in data:
            ncbi_lookup[item["AccessionVersion"]] = {
                "accession": item["AccessionVersion"],
                "id": item["Id"],
                "taxid": int(item["TaxId"]),
            }
        """
    return ncbi_lookup
