#! /usr/bin/env python
"""
Use a set of query reads to sweep out overlapping reads from another file.

% python sweepReads.py <ref_seqs> <reads> <nontarget_outfile> <target_outfile>

use "-" for <target_outfile> for stdout

Use '-h' for parameter help.
"""

import sys
import khmer
import numpy
import os.path
import time
import screed
from khmer.khmer_args import build_hashbits_args, report_on_config

def check_percentile(value):
    i = int(value)
    if not i in range(0, 101):
        raise argparse.ArgumentTypeError('%s must be between 0 to 100 inclusive' %value)
    return i

def main():
    parser = build_hashbits_args()
    parser.add_argument('ref_filename')
    parser.add_argument('read_filename')
    parser.add_argument('-P', '--kmer_identity',
                       default=50,
                       type=check_percentile,
                       help='kmer percent identity cutoff; '\
                            'percent of kmers in read covered by database')
    parser.add_argument('nontarget_output', 
                       help='nontarget reads output file')
    parser.add_argument('target_output',
                       help='target reads output file; use "-" for stdout')

    args = parser.parse_args()

    report_on_config(args, hashtype='hashbits')

    K = args.ksize
    HT_SIZE = args.min_tablesize
    N_HT = args.n_tables

    inp = args.ref_filename
    readsfile = args.read_filename
    perc = 100 - args.kmer_identity

    if args.target_output == '-':
        outfp = sys.stdout
    else:
        outfp = open('%s' %args.target_output, 'wb')

    nontarget_outfp = open('%s' %args.nontarget_output, 'wb')

    pt_file = '%s.pt' %(inp)

    if os.path.isfile(pt_file) and os.path.getsize(pt_file):
        # if hashbits table (.pt) exits
        ht = khmer.load_hashbits(pt_file)

    else:
        # load new hashbits table
        # create a hashbits data structure
        ht = khmer.new_hashbits(K, HT_SIZE, N_HT)

        # load contigs, connect into N partitions
        print >> sys.stderr, 'loading input reads from %s' %inp
        ht.consume_fasta(inp)

        print 'saving k-mer presence table in %s' %(pt_file)
        ht.save(pt_file)
 
    # Change 0.2 only if you really grok it.  HINT: You don't.
    fp_rate = khmer.calc_expected_collisions(ht)
    print >> sys.stderr, 'fp rate estimated to be %1.3f' % fp_rate

    if fp_rate > 0.20:
        print >>sys.stderr, "**"
        print >>sys.stderr, "** ERROR: the counting hash is too small for"
        print >>sys.stderr, "** this data set.  Increase hashsize/num ht."
        print >>sys.stderr, "**"
        print >>sys.stderr, "** Do not use these results!!"
        sys.exit(-1)

    print >> sys.stderr, 'starting sweep.'

    n = 0
    m = 0
    totalBp = 0
    start = time.time()
    for record in screed.open(readsfile):
        seqLen = len(record.sequence)
        totalBp += seqLen
        if seqLen < K:
            continue

        if n % 1000000 == 0:
            print >> sys.stderr, '... %d %d' %(n, m)

        count = ht.get_percentile_count(record.sequence, perc)

        if count:
            m += 1
            outfp.write('>%s\n%s\n' % (record.name, record.sequence))
        else:
            nontarget_outfp.write('>%s\n%s\n' % (record.name, record.sequence))
        n += 1

    end = time.time()
    duration = end - start
    print >> sys.stderr, 'Stats of the run:'
    print >> sys.stderr, '%.1f seqs per second' %(n*1.0/duration)
    print >> sys.stderr, '%.1f Mbps per second' %(totalBp*(1e-6)/duration)

if __name__ == '__main__':
    main()
