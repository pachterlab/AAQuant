import os
from math import log

def get_pvalue_weights(infile, outfile):

    if os.stat(infile).st_size == 0:
        open(outfile, 'a').close()
        return

    log_means = {}
    with open(infile, 'r') as fh:
        _, old_kmer, acc = fh.readline().rstrip('\n').split()
        acc = float(acc)
        N = 1
        for line in fh:
            _, kmer, count = line.rstrip('\n').split()
            try:
                count = float(count)
            except Exception as e:
                print(e)
                print(count)
                return
            if kmer != old_kmer:
                log_means[old_kmer] = log((acc/N) + 1)
                acc = 0.
                N = 0
                old_kmer = kmer
            if count > 0:
                acc += count
                N += 1

    weights = {kmer: log_m/sum(log_means.values()) for kmer, log_m in log_means.items()}

    with open(outfile, 'w') as fh:
        fh.write('\n'.join([f'{km}\t{w}' for km, w in weights.items()]))

if __name__ == '__main__':
    import sys
    get_pvalue_weights(sys.argv[1], sys.argv[2])
