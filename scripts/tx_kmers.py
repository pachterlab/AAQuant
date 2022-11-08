#!/usr/bin/env python3

import argparse
import re
from subprocess import check_output

K = 31

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
    # Who the heck puts semicolons in a semicolon-delimited file?!
    PATTERN = re.compile(r'''((}:[^;"']|"[^"]*"|'[^']*')+)''')
    rest = PATTERN.split(rest.rstrip('\n'))[1::2]
    try:
        return dict(map(lambda x: x.rstrip('\n').lstrip(' ').replace('"', '')\
                                                            .split(' ', 1), rest))
    except ValueError as e:
        print(rest)

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
    return scaffolds

def parse_gtf(gtf_path, fa_path, outfile, chrom, lb, ub):
    t_kmers = set()
    fa = parse_fa(fa_path)
    with open(gtf_path, 'r') as fh:
        for line in fh:
            # Skip comments
            if line.startswith('#'):
                continue

            data = line.split('\t')
            fields = {
                'chrom': data[0],
                'feature': data[2],
                'lb': int(data[3]) - 1,
                'ub': int(data[4]) - 1,
                'strand': data[6],
                'rest': parse_rest(' '.join(data[8]))
            }

            if fields['feature'] not in ['exon', 'UTR', 'stop_codon']:
                continue
            # Range selected
            if chrom is not None:
                if not fields['chrom'] == chrom:
                    continue
                # Make sure we're within [lb, ub]
                if fields['ub'] < lb or fields['lb'] > ub:
                    continue

            seq = fa[fields['chrom']][fields['lb']:fields['ub']]

            if fields['strand'] == '-':
                seq = reverse_complement(seq)
            t_kmers.update({seq[i:i+K] for i in range(len(seq)-K+1)})
    t_kmers = sorted(list(t_kmers))
    with open(outfile, 'w') as fh:
        fh.write('\n'.join(t_kmers))

def tx_kmers(gtf, fasta, gene, outfile):
    locus = check_output(['genecoords', '-w', gene]).decode('utf-8').rstrip('\n').split('\t')
    if locus[0] == '':
        return
    chrom, lb, ub = parse_range(locus[0])
    parse_gtf(gtf, fasta, outfile, chrom, lb, ub)

def tx_kmers_region(gtf, fasta, region, outfile):
    chrom, lb, ub = parse_range(region)
    parse_gtf(gtf, fasta, outfile, chrom, lb, ub)
