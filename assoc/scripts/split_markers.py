#!/usr/bin/env python3
import sys

def split_markers(in_path, out_path):
    with open(in_path, 'r') as fh:
        old_mk = ''
        for line in fh:
            pn, mk, a1, a2 = line.split()
            if old_mk != mk:
                try:
                    outfile.close()
                except NameError as e:
                    pass
                outfile = open(f'{out_path}/{mk}', 'w')
                old_mk = mk
            # TODO:
            # ====
            # Should probably replace the zeroes here with Minor Allele Frequency
            a1, a2 = max(0., float(a1)), max(0., float(a2))
            outfile.write(f'{pn}\t{a1+a2}\n')
        outfile.close()


if __name__ == '__main__':
    split_markers(sys.argv[1], sys.argv[2])
