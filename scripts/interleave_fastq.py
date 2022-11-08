#!/nfs_mount/bioinfo/apps-x86_64/python/3.7.3/bin/python3.7
import sys
label = 1
for line in sys.stdin:
    sys.stdout.write(f'>{label}\n')
    sys.stdout.write(line)
    label += 1
