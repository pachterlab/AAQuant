#!/usr/bin/env python3
import sys

def subset_counts(counts_path, subset):
    with open(counts_path, 'r') as fh:
        counts = [line.split() for line in fh]
    return ['\t'.join(c) for c in counts if c[1] in subset]

if __name__ == '__main__':
    with open(sys.argv[2], 'r') as fh:
        subset = sorted([line.split('\t')[0].rstrip('\n') for line in fh])
    print('\n'.join(subset_counts(sys.argv[1], subset)))
