#!/usr/bin/env python3
import argparse
import json
import os
import sqlite3
import subprocess
import sys

from collections import Counter
from functools import lru_cache

from sqlalchemy import create_engine, MetaData, Table, Integer, String, Column
from sqlalchemy.orm import mapper, sessionmaker
from sqlalchemy.pool import StaticPool


from taxcheck.lineage import LineageFactory, Lineage
from taxcheck.ncbi import ncbi_tax_lookup


class Gene:
    ...



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

@lru_cache(maxsize=10000)
def get_tax_annotation(db_session, gene_id):
    # print(f"Querying {gene_id} (type: {type(gene_id)}", file=sys.stderr)
    gene = db_session.query(Gene).filter(Gene.ACCESSION_VERSION == gene_id).one_or_none()
    return gene.TAXID if gene is not None else -1



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bamfile", type=str)
    # ap.add_argument("email", type=str)
    ap.add_argument("taxdb", type=str)
    ap.add_argument("--debug", action="store_true")
    ap.add_argument("--lineage-cutoff", type=float, default=0.75)
    ap.add_argument("--species-cutoff", type=float, default=0.99)
    # ap.add_argument("--ncbi-chunksize", type=int, default=400)
    # ap.add_argument("--ncbi_cache", type=str)
    args = ap.parse_args()

    if args.species_cutoff < args.lineage_cutoff:
        raise ValueError("Species cutoff ({args.species_cutoff}) needs to be greater than lineage cutoff ({args.lineage_cutoff}).")
    # if args.ncbi_chunksize < 10:
    #    raise ValueError("NCBI chunk size needs to be at least 10.")

    with open(args.taxdb, "rt") as db_in:
        db = {line.strip().split()[0]: tuple(line.strip().split()[1:]) for line in db_in}
    # source = sqlite3.connect(f"file:{args.taxdb}?mode=ro", uri=True)

    # # https://stackoverflow.com/questions/68286690/copy-an-sqlite-database-into-memory-using-sqlalchemy-for-testing-flask-app !!!
    # engine = create_engine("sqlite://", poolclass=StaticPool, connect_args={'check_same_thread': False},)
    # conn = engine.raw_connection().connection
        
    # source.backup(conn)
    
    # # engine = create_engine(f"sqlite:///{args.taxdb}")
    # metadata = MetaData(engine)

    # gene_table = Table(
    #     'NUCACC', metadata, 
    #     Column("ACCESSION_VERSION", String, primary_key=True),
    #     Column("TAXID", Integer),
    #     Column("GBID", Integer),
    # )
    # _ = mapper(Gene, gene_table)

    # Session = sessionmaker(bind=engine)
    # session = Session()

    cmd = ("samtools", "view", "-F", "0x4", args.bamfile)
    sam_proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    lfactory = LineageFactory()
    read2ref = {}
    # refs = set()
    ncbi_lookup = {}

    print("Parsing bam file...", file=sys.stderr)
    for nreads, (rname, ref, mapq, yp_tag, xa_tag) in enumerate(extract_yp_reads_from_sam(sam_proc.stdout), start=1):
        if nreads % 1000000 == 0:
            print(f"Processed {nreads} reads.", file=sys.stderr)
        
        note = ""
        lcount = Counter()
        refs = (ref, ) + (tuple(item[0] for item in xa_tag) if xa_tag else tuple())
        for gene_id in refs:
            acc =[x for x in gene_id.split("|") if x and x != "ref"]
            
            # lcount[ncbi_lookup.setdefault(gene_id, get_tax_annotation(session, acc[0]))] += 1
            lcount[ncbi_lookup.setdefault(gene_id, db.get(acc[0]))] += 1

        lineages = {
            taxid: {
                "lineage": lfactory.generate_lineage(taxid),
                "count": count,
            }
            for taxid, count in lcount.items()
        }

        if len(lineages) == 1:
            consensus_lineage = tuple(lineages.values())[0]["lineage"]
            consensus_level = "species"
            consensus_id, consensus_name = consensus_lineage.levels[-1].values()
        else:

            # iterate from species -> kingdom
            for level in range(Lineage.TAXLEVELS["species"][0], -1, -1):
                tax_counter = Counter()
                # for each level count the taxids (multiplied by number of alignments)
                for lineage_data in lineages.values():
                    lineage, count = lineage_data["lineage"], lineage_data["count"]
                    level_id = lineage.levels[level]["id"] if lineage.levels[level] is not None else -1
                    tax_counter.update((level_id,) * count)
                if args.debug:
                    print("LEVEL", level, tax_counter, file=sys.stderr)
                # then check if there's a consensus (based on fractional representation by the alignments)
                if tax_counter:
                    top_taxid, top_count = tax_counter.most_common()[0]
                    cutoff = args.species_cutoff if level == Lineage.TAXLEVELS["species"][0] else args.lineage_cutoff
                    if top_taxid == -1:
                        if not note:
                            note = f"FIRST_LEVEL_UNKNOWN={list(Lineage.TAXLEVELS)[level]}"
                    elif top_count / sum(tax_counter.values()) > cutoff:
                        consensus_lineage = lfactory.generate_lineage(top_taxid)
                        consensus_level = list(Lineage.TAXLEVELS)[level]
                        consensus_id, consensus_name = consensus_lineage.levels[level].values()
                        break


        # lcount[ncbi_lookup.get(ref, {}).get("taxid", -1)] += 1
        
        #     # count how often a read aligns against a refseq from a specific taxonomy
        #     lcount = Counter()
        #     for ref in rrefs:
        #         if args.debug:
        #             print(ref, ncbi_lookup.get(ref), file=sys.stderr)
        #         lcount[ncbi_lookup.get(ref, {}).get("taxid", -1)] += 1

        #     if args.debug:
        #         print("LCOUNT", lcount, file=sys.stderr)
        #         print("LCOUNT HAS -1", -1 in lcount, file=sys.stderr)

        # generate lineages from each reference and store the alignment counts along with it
        # lineages2 = {
        #     taxid: {
        #         "lineage": lfactory.generate_lineage(taxid),
        #         "count": lcount[taxid],
        #     }
        #     for taxid in yp_tag  # aln_data["yp"]
        # }

        # if args.debug:
        #     print(*lineages2.values(), sep="\n", file=sys.stderr)
        #     # print(*(f"{l.get_string()}: {c}" for l, c in lineages2.values()), sep="\n", file=sys.stderr)

        # determine a consensus line based on fraction of alignments
        # if len(lineages2) == 1:
        #     consensus_lineage = tuple(lineages2.values())[0]["lineage"]
        #     consensus_level = "species"
        #     consensus_id, consensus_name = consensus_lineage.levels[-1].values()
        # else:
        #     # iterate from species -> kingdom
        #     for level in range(Lineage.TAXLEVELS["species"][0], -1, -1):
        #         tax_counter = Counter()
        #         # for each level count the taxids (multiplied by number of alignments)
        #         for lineage_data in lineages2.values():
        #             lineage, count = lineage_data["lineage"], lineage_data["count"]
        #             level_id = lineage.levels[level]["id"] if lineage.levels[level] is not None else -1
        #             tax_counter.update((level_id,) * count)
        #         if args.debug:
        #             print("LEVEL", level, tax_counter, file=sys.stderr)
        #         # then check if there's a consensus (based on fractional representation by the alignments)
        #         if tax_counter:
        #             top_taxid, top_count = tax_counter.most_common()[0]
        #             cutoff = args.species_cutoff if level == Lineage.TAXLEVELS["species"][0] else args.lineage_cutoff
        #             if top_taxid == -1:
        #                 if not note:
        #                     note = f"FIRST_LEVEL_UNKNOWN={list(Lineage.TAXLEVELS)[level]}"
        #             elif top_count / sum(tax_counter.values()) > cutoff:
        #                 consensus_lineage = lfactory.generate_lineage(top_taxid)
        #                 consensus_level = list(Lineage.TAXLEVELS)[level]
        #                 consensus_id, consensus_name = consensus_lineage.levels[level].values()
        #                 break

        print(
            rname,
            len(lineages),
            # len(aln_data["xa"]) + 1,
            len(xa_tag) + 1,
            consensus_level,
            consensus_id,
            consensus_name,
            consensus_lineage.get_string(),
            consensus_lineage.get_string(show_names=False),
            note,
            sep="\t",
        )


if __name__ == "__main__":
    main()
