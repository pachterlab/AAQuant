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
    # Who the heck puts semicolons in a semicolon-delimited file?!
    PATTERN = re.compile(r'''((}:[^;"']|"[^"]*"|'[^']*')+)''')
    rest = PATTERN.split(rest.rstrip('\n'))[1::2]
    try:
        return dict(map(lambda x: x.rstrip('\n')\
                                   .lstrip(' ')/
                                   .replace('"', '')\
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

def parse_gtf(gtf_path, fa_path, chrom, lb, ub):
    exons = []
    fa = parse_fa(fa_path)
    with open(gtf_path, 'r') as fh:
        for line in fh:
            # Skip comments
            if line.startswith('#'):
                continue

            data = line.split('\t')
            fields = {
                'chrom'   : data[0],
                'feature' : data[2],
                'lb'      : int(data[3]) - 1,
                'ub'      : int(data[4]) - 1,
                'strand'  : data[6],
                'rest'    : parse_rest(' '.join(data[8]))
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
            exons.append(seq)
    return exons

def get_exons(gtf, fasta, gene):
    locus = check_output(['genecoords', '-w', gene]).decode('utf-8').rstrip('\n').split('\t')
    if locus[0] == '':
        return
    chrom, lb, ub = parse_range(locus[0])
    return parse_gtf(gtf, fasta, chrom, lb, ub)

def parse_k2u(path):
    k2u = {}
    with open(path, 'r') as fh:
        for line in fh:
            k, u = line.split()
            k2u[k] = u.rstrip('\n')
    return k2u

def detect_splice(kmer, exons):
    km = kmer[:16]
    er = kmer[16:]

    hits = []

    km_loc = -1
    er_loc = -1

    for exon, seq in enumerate(exons):
        # First pass to look for one(plus) half of kmer

        # TODO:
        # Look for reverse complement as well!
        km_loc = seq.find(km)
        er_loc = seq.find(er)

        if km_loc == er_loc == -1 or er_loc == (km_loc+len(km)):
            # a) Neither km nor er in current exon (boo!)
            # b) Exon contains km and er as a whole (boring)
            continue

        # Expand the exact matching of the subsequence that was found
        km_exp = 0
        if km_loc > -1:

            for i, e in enumerate(er):
                if len(seq) > km_loc + len(km) + i and seq[km_loc + len(km) + i] == e:
                    km_exp += 1
                else:
                    break

            if len(km) + km_exp == len(kmer):
                continue

            hits.append({
                'km'      : km + er[:km_exp],
                'km_exon' : exon,
                'km_loc'  : km_loc,
                'er'      : er[km_exp:],
                'er_exon' : -1,
                'er_loc'  : -1
            })

        if er_loc > -1:

            if km_exp == 0 and km_loc > -1:
                # km and er are split exactly [:16] and [16:], respectively
                hits[-1]['er_exon'] = exon
                hits[-1]['er_loc']  = er_loc

            else:
                # Check whether a larger er exists in exon
                er_exp = 0
                for i in range(1, len(km) + 1):
                    if (er_loc - i) > -1 and seq[er_loc - i] == km[-i]:
                        er_exp += 1
                    else:
                        break

                if len(er) + er_exp == len(kmer):
                    continue

                hits.append({
                    'km'      : km[:-er_exp] if er_exp > 0 else km,
                    'km_exon' : -1,
                    'km_loc'  : -1,
                    'er'      : km[-er_exp:] + er if er_exp > 0 else er,
                    'er_exon' : exon,
                    'er_loc'  : er_loc - er_exp
                })

    for exon, seq in enumerate(exons):
        # Second pass to look for other half of kmer

        # TODO:
        # ACTUALLY IMPLEMENT THIS PROPERLY!
        # I think check?
        for hit in hits:
            # Gimme that O(m * n^2)
            if hit['km_loc'] < 0:

                km_len = len(hit['km'])
                if km_len < 5:
                    continue

                km_loc = seq.find(hit['km'])
                if km_loc < 0:
                    continue

                os = max(len(seq), km_loc+km_len+3)
                if len(seq) - km_loc - km_len <= 0:
                    hit['km_exon'] = exon
                    hit['km_loc']  = km_loc

                elif seq[km_loc + km_len:os] != hit['er'][:os-km_loc-km_len]:
                    # len(seq) > km_loc + km_len
                    hit['km_exon'] = exon
                    hit['km_loc']  = km_loc

            elif hit['er_loc'] < 0:

                if len(hit['er']) < 5:
                    continue

                er_loc = seq.find(hit['er'])
                if er_loc < 0:
                    continue

                os = min(er_loc, 3)

                if er_loc == 0 or seq[er_loc - os:er_loc] != hit['km'][-os:]:
                    hit['km_exon'] = exon
                    hit['km_loc']  = er_loc


    return hits

def process_kmers():

    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gtf', type=str, help='Path to GTF file.')
    parser.add_argument('-f', '--fasta', type=str, help='Path to FASTA file.')
    parser.add_argument('--gene', type=str, help='Name of gene.')
    parser.add_argument('-a', '--assoc', type=str, help='Path to associations.')
    parser.add_argument('--k2u', type=str, help='kmer to unitig mapping.')
    parser.add_argument('--exclude-retention', dest='exclude', action='store_true')
    args = parser.parse_args()

    exons = get_exons(args.gtf, args.fasta, args.gene)

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

            hits = detect_splice(unitig, exons) + detect_splice(reverse_complement(unitig), exons)

            for hit in hits:
                #if not args.exclude or (hit['km_exon'] > -1 and hit['er_exon'] > -1):
                #                                            XOR
                if not args.exclude or (hit['km_exon'] == -1) != (hit['er_exon'] == -1):
                    print(unitig, '\t'.join([str(h) for h in hit.values()]), '\t'.join(rest), sep='\t')

if __name__ == '__main__':
    process_kmers()
