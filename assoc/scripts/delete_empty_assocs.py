#!/usr/bin/env python3

import os, sys
from glob import iglob

def has_one_line(f):
    with open(f, 'r') as fh:
        return fh.readline() == 'kmer\tcoef\tstderr\tdf\n' and fh.readline() == ''


if __name__ == '__main__':
    root = '/odinn/tmp/kristjanh/output_11062020/*/assoc/*'
    if len(sys.argv) > 1:
        root = sys.argv[1]
    files = (f for f in iglob(root, recursive=True) if os.path.isfile(f))
    deleted = 0
    for f in files:
        if has_one_line(f):
            print(f'rm {f}')
            os.remove(f)
            deleted += 1
    print(f'removed {deleted} files.')

