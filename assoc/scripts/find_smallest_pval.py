#!/usr/bin/env python3
import os, sys
from glob import iglob


if __name__ == '__main__':
    root = '/odinn/tmp/kristjanh/output_11062020/*/assoc/*'
    if len(sys.argv) > 1:
        root = sys.argv[1]
    files = (f for f in iglob(root, recursive=True) if os.path.isfile(f))

    assocs = {}
    for f in files:
        with open(f, 'r') as fh:
            fh.readline()
            for line in fh:
                kmer, _, _, _, p = line.split()
                p = float(p)
                if f in assocs and assocs[f]['p'] > p:
                    assocs[f]['p'] = p
                    assocs[f]['kmer'] = kmer
                elif f not in assocs:
                    assocs[f] = {'p': p, 'kmer': kmer, 'path': f}

    asc = sorted(assocs.values(), key=lambda e: e['p'])

    for a in asc[:40]:
        print(f'p: {a["p"]}, kmer: {a["kmer"]}, in file {a["path"]}')
