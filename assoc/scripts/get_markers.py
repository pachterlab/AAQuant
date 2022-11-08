#!/usr/bin/env python3
from subprocess import run, check_output as co
from os import makedirs

if __name__ == '__main__':
    genes = [line.rstrip('\n') for line in open('/nfs/odinn/users/kristjanh/work/kmer_association/dbg/prune_pipeline/gene_lists/100genes.txt') if not line.startswith('#')]
    for gene in genes:
        rg, _ = co(['genecoords', '-w', gene]).decode('utf-8').rstrip('\n').split()
        chrom, iv = rg.split(':')
        iv = '-'.join(str(int(t[0])+t[1]) for t in zip(iv.split('-'), [-40000, 40000]))
        makedirs(f'../markers/{gene}', exist_ok=True)
        run(' '.join(['chitools',
               'stats',
               '-p',
               '/nfs/odinn/users/kristjanh/work/kmer_association/dbg/prune_pipeline/pn_lists/blood_bamfiles.filtered.pns',
               '--interval',
               iv,
               f'/nfs/odinn/data/results/imputations/aster/imputedseq/{chrom}.chi',
               '|',
               'awk',
               '\'{if ($6>0.01) print $2}\'',
               '|',
               'tail',
               '-n',
               '+2',
               '>',
               f'../markers/{gene}/marker_names_40kb_freq_gt.01']), shell=True)
        run(' '.join(['chitools',
               'view',
               '-p',
               '/nfs/odinn/users/kristjanh/work/kmer_association/dbg/prune_pipeline/pn_lists/blood_bamfiles.filtered.pns',
               '-m',
               f'../markers/{gene}/marker_names_40kb_freq_gt.01',
               f'/nfs/odinn/data/results/imputations/aster/imputedseq/{chrom}.chi',
               '>',
               f'../markers/{gene}/markers_40kb_freq_gt.01']), shell=True)
