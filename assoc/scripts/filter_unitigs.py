#!/usr/bin/env python3

import argparse
import re
from subprocess import check_output

def reverse_complement(string):
    comp = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        'M': 'K',
        'K': 'M',
        'Y': 'R',
        'R': 'Y',
        'S': 'S',
        'W': 'W',
        'N': 'N'
    }
    return ''.join([comp[c] for c in string[::-1]])

def parse_range(rstring):
    chrom, rest = rstring.split(':')
    lb, ub = (int(i) for i in rest.split('-'))
    return chrom, lb, ub

def parse_rest(rest):
    return dict([r.split(' ') for r in rest.rstrip(';\n').replace('"', '').split('; ')])

def parse_fa(fa_path):
    scaffolds = {}
    scaffold = ''
    split_lines = []
    with open(fa_path, 'r') as fh:
        for line in fh:
            if line.startswith('>'):
                if scaffold != '':
                    scaffolds[scaffold] = ''.join(split_lines)
                scaffold = line.split()[0][1:]
                split_lines = []
            else:
                split_lines.append(line.rstrip('\n'))
        scaffolds[scaffold] = ''.join(split_lines)
    return scaffolds

def parse_gtf(gtf_path, fa_path, gene):
    transcripts = []
    fa = parse_fa(fa_path)
    with open(gtf_path, 'r') as fh:
        for line in fh:
            # Skip comments
            if line.startswith('#'):
                continue

            data = line.split('\t')
            f = {
                'chrom'   : data[0],
                'feature' : data[2],
                'lb'      : int(data[3]) - 1,
                'ub'      : int(data[4]) - 1,
                'strand'  : data[6],
                'rest'    : parse_rest(''.join(data[8]))
            }

            if f['feature'] != 'transcript' or f['rest']['gene_name'] != gene or f['rest']['gene_biotype'] != 'protein_coding':
                continue

            """
            # Range selected
            if chrom is not None:
                if not f['chrom'] == chrom:
                    continue
                # Make sure we're within [lb, ub]
                if f['ub'] < lb or f['lb'] > ub:
                    continue
            """

            # In out .fa file, reverse strand transcripts are already reverse complemented
            seq = fa[f"{f['rest']['transcript_id']}.{f['rest']['transcript_version']}"]

            transcripts.append(seq)
    return transcripts

"""
def get_exons(gtf, fasta, gene):
    locus = check_output(['genecoords', '-w', gene]).decode('utf-8').rstrip('\n').split('\t')
    if locus[0] == '':
        return
    chrom, lb, ub = parse_range(locus[0])
    return parse_gtf(gtf, fasta, chrom, lb, ub)
"""

def parse_k2u(path):
    k2u = {}
    with open(path, 'r') as fh:
        for line in fh:
            k, u = line.split()
            k2u[k] = u.rstrip('\n')
    return k2u

def filter_unitigs(kmer, exons):
    for e in exons:
        if kmer in e:
            return False
    return True


def process_kmers():

    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gtf', type=str, help='Path to GTF file.')
    parser.add_argument('-f', '--fasta', type=str, help='Path to FASTA file.')
    parser.add_argument('--gene', type=str, help='Name of gene.')
    parser.add_argument('-a', '--assoc', type=str, help='Path to associations.')
    parser.add_argument('--k2u', type=str, help='kmer to unitig mapping.')
    parser.add_argument('--exclude-retention', dest='exclude', action='store_true')
    args = parser.parse_args()

    transcripts = parse_gtf(args.gtf, args.fasta, args.gene)

    if len(args.k2u) > 0:
        k2u = parse_k2u(args.k2u)

    with open(args.assoc, 'r') as fh:

        fh.readline()
        for line in fh:

            kmer, *rest = line.split()

            if len(args.k2u) > 0:
                unitig = k2u[kmer]
            else:
                unitig = kmer

            if filter_unitigs(unitig, transcripts) and filter_unitigs(reverse_complement(unitig), transcripts):
                print(line.rstrip('\n'))

if __name__ == '__main__':
    process_kmers()
