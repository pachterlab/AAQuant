import sys
from subprocess import check_output

def get_locus(gene):
    try:
        locus, _ = check_output(['genecoords', '-w', gene]).decode('utf-8').rstrip('\n').split('\t')
        return locus
    except ValueError as e:
        return None

def print_tab(*args):
    print('\t'.join(args))

if __name__ == '__main__':
    events = []
    with open(sys.argv[1], 'r') as fh:
        for line in fh:
            gene = line.rstrip('\n')
            locus = get_locus(gene)
            chrom, rest = locus.split(':')
            lb, ub = rest.split('-')
            events.append({'gene': gene, 'chr': chrom, 'lb': lb, 'ub': ub})

    with open(sys.argv[2], 'r') as fh:
        fh.readline()
        for line in fh:
            chrom, lb, _, sid, _, svlen, *rest = line.split()
            for e in events:
                if chrom == e['chr']:
                    if e['lb'] <= lb <= e['ub'] or e['lb'] <= lb + svlen <= e['ub']:
                        print_tab(e['gene'], sid, lb, svlen)
