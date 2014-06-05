#! usr/bin/python
# chang seq names in fasta file
# by gjr; Mar 21, 11

import sys,os
import screed

def main():
    # just index seqs
    if len(sys.argv) != 4:
        print >> sys.stderr, \
                     'Usage: python %s <seqfile><out.fasta><tag>' \
                     %(os.path.basename(sys.argv[0]))
        sys.exit(1)

    tag = sys.argv[3]
    fw = open(sys.argv[2], 'w')
    for n, record in enumerate(screed.open(sys.argv[1])):
        name = record['name']
        seq = record['sequence']
        name = name.strip()
        _lis = name.split()
        if len(_lis) != 1:
            print >> sys.stderr, 'name with space dectected %s' %(name)
            print >> sys.stderr, "replacing ' ' with '_'"
            name = '_'.join(_lis)

        fw.write('>%s%d__%s\n%s\n' %(tag, n, tag, seq))

if __name__ == '__main__':
    main()
