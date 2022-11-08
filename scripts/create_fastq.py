from subprocess import check_output

def create_fastq(gene, samples, outfile):
    locus, _  = check_output(['genecoords', '-w', gene]).decode('utf-8').rstrip('\n').split('\t')
    with open(outfile, 'w') as of:
        for pn, bam in samples.items():
            out = check_output(['samtools', 'view', bam, locus])
            try:
                reads = [o.split()[9] for o in out.decode().rstrip('\n').split('\n')]
                for i, r in enumerate(reads):
                    of.write(f'>{pn}.{i}\n')
                    of.write(f'{r}\n')
            except IndexError as e:
                print(f'{pn} does not express {gene}, {locus}')
