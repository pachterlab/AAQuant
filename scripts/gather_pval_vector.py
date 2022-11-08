#!/usr/bin/env python3
import os
import matplotlib.pyplot as plt

if __name__ == '__main__':
    import sys

    path = sys.argv[1]
    line = int(sys.argv[2])
    ticks = []
    pvals = []
    for f in os.listdir(path):
        with open(os.path.join(path, f)) as fh:
            lines = fh.readlines()
            ticks.append(float(f.split('_')[-1]))
            pvals.append(float(lines[line].rstrip('\n').split()[-1]))

    print(ticks)
    print(pvals)
    pvals = [x for y, x in sorted(zip(ticks, pvals))]
    ticks.sort()
    plt.plot(ticks, pvals)
    plt.yscale('log')
    plt.plot(ticks, [0.05 for _ in ticks], 'g--')
    # plt.text(ticks[-5], 100, 'p=0.05', color='g')
    plt.text(ticks[-5], 0.06, 'p=0.05', color='g')
    plt.plot(ticks, [0.05/20000 for _ in ticks], 'r--')
    # plt.text(ticks[-5], 1e-15, 'p=2.5e-6', color='r')
    plt.text(ticks[0], 3e-6, 'p=2.5e-6', color='r')
    plt.xlabel('True effect size')
    plt.ylabel('log(p)')
    plt.title('p-value of structurally significant contig')
    plt.show()
