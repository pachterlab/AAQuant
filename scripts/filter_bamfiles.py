#WHITELIST_PATH = '/odinn/data/dataprocessing/rnaexpr-adipose/bird/rnacoverage/filtered.coverage.SETID.txt'
WHITELIST_PATH = '/odinn/data/dataprocessing/rnaexpr-blood/bird/rnacoverage/filtered.coverage.SETID.txt'
#BAMFILES = 'adipose_bamfiles.txt'
BAMFILES = 'blood_bamfiles.txt'

if __name__ == '__main__':
    whitelist = [w.rstrip('\n') for w in open(WHITELIST_PATH, 'r')]
    with open(BAMFILES, 'r') as fh:
        for line in fh:
            setid, path = line.rstrip('\n').split()
            if setid in whitelist:
                print(line)
