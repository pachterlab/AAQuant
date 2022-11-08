import argparse
import re
from operator import itemgetter

def parse_fasta(path, s_filter=None):
    scaffolds = {}
    sequences = []
    current = ''
    with open(path, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if len(current) > 0 and (s_filter is None or current in s_filter):
                    scaffolds[current] = ''.join(sequences)
                    if s_filter is not None:
                        break
                sequences = []
                current = line.split()[0][1:]
            elif current is None or current in s_filter:
                sequences.append(line)
        scaffolds[current] = ''.join(sequences)
    return scaffolds

def merge_intervals(intervals):
    # Sort by lower bound
    sorted_intervals = sorted(intervals, key=itemgetter(0))
    # No intervals
    if not sorted_intervals:
        return
    low, high = sorted_intervals[0]
    for iv in sorted_intervals[1:]:
        if iv[0] <= high:
            high = max(high, iv[1])
        else:
            yield low, high
            low, high = iv
    yield low, high

def parse_rest(rest):
    # Parses the 'Attributes' field of the gtf
    # Who the heck puts semicolons in a semicolon-delimited file?!
    PATTERN = re.compile(r'''((?:[^;"']|"[^"]*"|'[^']*')+)''')
    return dict(map(lambda x: x.rstrip('\n')\
                               .lstrip(' ')\
                               .replace('"', '')\
                               .split(' ', 1), rest))

def parse_gtf(gtf_path, gene):
    trs = {}
    with open(gtf_path, 'r') as fh:
        for line in fh:
            if line.startswith('#'):
                continue

            data = line.split('\t')
            fields = {
                'scaffold': data[0],
                'feature' : data[2],
                'start'   : int(data[3]) - 1, # Scaffold is 0-indexed
                'end'     : int(data[4]) - 1,
                'strand'  : data[6],
                'rest'    : parse_rest(data[8])
            }
            gene_name = fields['rest']['gene_name']
            if gene_name != gene:
                continue

            if fields['feature'] != 'gene':
                tr = fields['rest']['transcript_id']
                if tr not in trs:
                    trs[tr] = {
                        'exon': [],
                        'stop_codon': [],
                        'utr': [],
                        'transcript': []
                    }
                if fields['feature'] in ['UTR', 'utr', 'five_prime_utr', 'three_prime_utr']:
                    trs[tr]['utr'].append(fields)
                else:
                    trs[tr][fields['feature']].append(fields)
    return trs

def generate_ref_filter(ref_fa, ref_gtf, transcriptome_fa, gene, output):
    trs = parse_gtf(ref_gtf, gene)

    # Find the scaffold we need to subset in the fasta file
    scaffold = None
    for features in trs.values():
        for feature in features.values():
            if feature:
                scaffold = feature[0]['scaffold']
                break
        if scaffold is not None:
            break

    fa = parse_fasta(ref_fa, [scaffold])

    sequences = {}
    features = []
    for trans, comp in trs.items():
        for typ, fts in comp.items():
            #ivs = merge_intervals(((f['start'], f['end']) for f in features))
            ivs = ((f['start'], f['end']) for f in fts)
            
            for i, f in enumerate(fts):
                version = f['rest']['transcript_version']
                name = f'{trans}.{version} {typ} {str(i)}'
                sequences[name] = fa[scaffold][f['start']:f['end']]
                features.append(f'{trans}.{version}')

    features = list(set(features))

    sequences = sequences.update(parse_fasta(transcriptome_fa, features))
    return sequences

def write_fasta(sequences, output):
    with open(output, 'w') as fh:
        for feature, seq in sequences.items():
            fh.write(f'>{feature}')
            fh.write(seq)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate ref filter from fasta and gtf files.')
    parser.add_argument('--ref-fa', type=str, help='Genome fasta file.')
    parser.add_argument('--ref-gtf', type=str, help='Genome gtf file.')
    parser.add_argument('--transcriptome-fa', type=str, help='Transcriptome fasta file.')
    parser.add_argument('--gene', type=str, help='Gene name.')
    parser.add_argument('--output', type=str, help='The file to which the ref filter should be written.')
    args = parser.parse_args()
    write_fasta(generate_ref_filter(args.ref_fa,
                                    args.ref_gtf,
                                    args.transcriptome_fa,
                                    args.gene,
                                    args.output))
