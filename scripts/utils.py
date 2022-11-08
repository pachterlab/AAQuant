import sys
from subprocess import check_output

def get_locus(gene):
    try:
        locus, _ = check_output(['genecoords', '-w', gene]).decode('utf-8').rstrip('\n').split('\t')
        return locus
    except ValueError as e:
        return None

if __name__ == '__main__':
    for line in sys.stdin:
        if get_locus(line.rstrip('\n')) is not None:
            sys.stdout.write(line)
