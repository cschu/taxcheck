#!/usr/bin/env python3
import argparse
import gzip
import subprocess
import sys




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




def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bamfile", type=str)
    args = ap.parse_args()

    read2ref = {}

    cmd = ("samtools", "view", "-F", "0xf04", args.bamfile)
    sam_proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    for aln in get_lines_from_chunks(sam_proc.stdout):
        flag = int(aln[1])
        if flag & 0x1 and not (aln[2] == aln[6] or flag & 0x8):
            continue

        tags = dict(item.split(":")[0::2] for item in aln[11:])
        if not tags.get("XA"):
            read2ref.setdefault(ref[1], set()).add(ref[0])              


            
            
    for ref in read2ref:
        print(ref)




if __name__ == "__main__":
    main()
