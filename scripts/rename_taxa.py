#! usr/bin/python
# chang seq names in fasta file
# by gjr; Mar 21, 11

import sys,os

def main():
    if len(sys.argv) != 4:
        print >> sys.stderr, \
                     'Usage: python %s <file.taxa><out.taxa><tag>' \
                     %(os.path.basename(sys.argv[0]))
        sys.exit(1)

    tag = sys.argv[3]
    with open(sys.argv[2], 'w') as fw:
        with open(sys.argv[1]) as fp:
            for line in fp: 
                name, taxa = line.rstrip().split('\t',1)
                name = name.strip()
                _lis = name.split()
                if len(_lis) != 1:
                    print >> sys.stderr, 'name with space dectected %s' %(name)
                    print >> sys.stderr, "replacing ' ' with '_'"
                    name = '_'.join(_lis)
                fw.write('%s__%s\t%s\n' %(name, tag, taxa))

if __name__ == '__main__':
    main()
