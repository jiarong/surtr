#! /usr/bin/env python
# calculate the edit_dist(NM)/alignmentLen(M and I in CIGAR)
# by gjr; July 25, 13


"""
M alignment match (can be a sequence match or mismatch)
I insertion to the reference
D deletion from the reference
N skipped region from the reference
S soft clipping (clipped sequences present in SEQ)
H hard clipping (clipped sequences NOT present in SEQ)
P padding (silent deletion from padded reference)
= sequence match
X sequence mismatch

CIGAR operation function from sqt by Marcel Martin.

There are two ways to represent a CIGAR string:
- as a string, such as "17M1D5M4S"
- as a list of (length, operator) pairs, as used by pysam:
[ (17, 0), (2, 1), (5, 0), (4, 0) ]

The naming convention in this module uses cigar and cigar_string to
distinguish both types.

The mapping of CIGAR operator to numbers is:
MIDNSHP => 0123456
"""
import sys
from itertools import repeat, chain
import fileinput
import argparse

# constants
M = 0
I = 1
D = 2
N = 3
S = 4
H = 5
P = 6
# assuming "X" and "=" are not used

# use this as a sequence to map an encoded operation to the appropriate
# character
DECODE = 'MIDNSHP'

# this dictionary maps operations to their integer encodings
_ENCODE = dict( (c,i) for (i, c) in enumerate(DECODE) )


def parse(cigar_string):
    """
    Parse CIGAR string and return a list of (length, operator) pairs.

    >>> parse("3S17M8D4M9I3H")
    [(4, 3), (0, 17), (2, 8), (0, 4), (1, 9), (5, 3)]
    """
    cigar = []
    n = ''  # this is a string, to which digits are appended
    for c in cigar_string:
        if c.isdigit():
            n += c
        elif c in _ENCODE:
            if n == '':
                raise ValueError("end of CIGAR string reached, but an operator was expected")
            cigar.append( (_ENCODE[c], int(n)) )
            n = ''
    return cigar


def as_string(cigar):
    """
    Convert CIGAR, given as list of (length, operator) pairs, to a string.

    >>> as_string([0,15), (2,1), (0,36)]
    M15D1M36
    """
    return ''.join('{}{}'.format(l, DECODE[op]) for op, l in cigar)


def concat(left, right):
    """
    Concatenate two CIGARs given as list of (length, operator) pairs.

    >>> concat(parse_cigar("1M"), parse_cigar("3M"))
    [(0, 4)]
    """
    if left_cigar and right:
        left_last = left[-1]
        right_first = right[0]
        # same operation?
        if left_last[0] == right_first[0]:
            right[0] = ( right[0][0], right[0][1] + left[-1][1] )
            left = left[:-1]
    return left + right


def aligned_length(cigar):
    """
    Return aligned length of a CIGAR given as list of (length, operator) pairs.
    """
    length = 0
    for op, l in cigar:
        if op in [0, 1]: # MI
            length += l
        elif op == 2: # D
            pass
        else:
            raise ValueError("CIGAR operation %s not supported" % cigar)
    return length


def ops(cigar):
    """
    Yield all operations (as numbers, not characters) one by one.

    >>> list(ops(parse("3S2I3M")))
    [4, 4, 4, 1, 1, 0, 0, 0]
    """
    return chain.from_iterable(repeat(op, l) for (op, l) in cigar)


def decoded_ops(cigar):
    """
    Yield all operations (as characters) one by one.

    >>> ''.join(ops(parse("3S2I3M")))
    "SSSIIMMM"
    """
    return chain.from_iterable(repeat(DECODE[op], l) for (op, l) in cigar)


def _assert_at_end(i):
    """Assert that the iterator i is at its end"""
    if __debug__:
        try:
            next(i)
            assert False
        except StopIteration:
            pass


def alignment_iter(read, ref, cigar, gap='-'):
    """
    Yield triples (read_char, reference_char, cigar_char) that
    fully describe the alignment betwen read and ref according to cigar.

    If the cigar operation is a 'M', the cigar_char is set to either
    '=' or 'X' depending on whether read_char matches reference_char
    or not.

    At gaps in the alignment, either read_char or reference_char are
    set to the given gap character.

    read -- an iterable representing the read
    ref -- an iterable representing the reference sequence
    cigar -- a list of (operator, length) pairs
    """
    i = iter(read)
    j = iter(ref)
    for op in decoded_ops(cigar):
        if op == 'M':
            ci = chr(next(i))
            cj = chr(next(j))
            yield (ci, cj, '=' if ci == cj else 'X')
        elif op == 'I':
            yield (chr(next(i)), gap, 'I')
        elif op == 'D':
            yield (gap, chr(next(j)), 'D')
        else:
            raise ValueError("CIGAR operator {} not supported".format(op))
    _assert_at_end(i)
    _assert_at_end(j)


def print_alignment(read, ref, cigar, file=sys.stdout):
    """
    Print an alignment between read and ref according to a CIGAR.
    This uses the alignment_iter() function from above.

    cigar -- a list of (operator, length) pairs
    """
    row1 = ''
    row2 = ''
    align = ''
    for read_char, reference_char, op in alignment_iter(read, ref, cigar):
        row1 += read_char
        align += op
        row2 += reference_char
    #print(row1, align, row2, sep='\n', file=file)


def unclipped_region(cigar):
    """
    Return tuple (cigar, start, stop), where cigar is the given cigar without soft clipping
    and (start, stop) is the interval in which the read is *not* soft-clipped.
    """
    if cigar[0][0] == S:
        start = cigar[0][1]
        cigar = cigar[1:]
    else:
        start = 0
    if cigar[-1][0] == S:
        stop = -cigar[-1][1]
        cigar = cigar[:-1]
    else:
        stop = None
    return (cigar, start, stop)

def check_positive(value):
    i = int(value)
    if i < 0:
        raise argparse.ArgumentTypeError('%s must be non-negative value' %value)
    return i

def main():
    # usage: %python <thisFile> -c cov -i id -l length -o outdir -f <samfileFormMEM>
    # replace <listOfReadNameFile> with - if pipe
    # e.g. bwa mem ref.fa reads.fa |python [options] <thisFile> - outdir
    # output should be prefix.rna, prefix.nonrna, prefix.csv

    parser = argparse.ArgumentParser(description='Parser to filter bwa output'\
                                                    'based on coverage' \
                                                    ' and identity') 
    
    parser.add_argument('-c', '--coverage', 
                       default=0,
                       type=check_positive,        
                       help='coverage percentage cutoff'\
                               '(100*alignment/read length);'\
                               ' Coverage of alignment on a read.')
    parser.add_argument('-i', '--align_identity', 
                       default=0,
                       type=check_positive,
                       help='alignment percentage identity cutoff '\
                               '(100*match/alignment length);'\
                               ' match percentage on the alignment.')
    parser.add_argument('-I', '--read_identity',
                       default=0,
                       type=check_positive,
                       help='read percentage identity cutoff '\
                               '(100*alignment/read length);'\
                               ' match percentage on the read.')
    parser.add_argument('-m', '--match_number',
                       default=0,
                       type=check_positive,
                       help='number of matches required.')
    parser.add_argument('-o', '--outdir',
                       default='./',
                       help='output directory')
    parser.add_argument('-p', '--prefix',
                       default='temp',
                       help='prefix for output filenames')

    parser.add_argument('inputfile',
                       help='input samfile; use "-" for stdin')

    args = parser.parse_args()


    with open('%s/%s.nonrrna' %(args.outdir, args.prefix), 'wb') as fw_nonrna,\
         open('%s/%s.temp' %(args.outdir, args.prefix), 'wb') as fw_temp:

        d = {}  # dict to collect info
        d_ref_length = {}
         
        nm_tag_pos = None
        for n, line in enumerate(fileinput.input(args.inputfile)):
            # for bwa-mem:
            # qname, flag, rname, pos, mapQ, cigar, rnext, pnext, 
            #    tlen, seq, qual, NM, others
            lis = line.rstrip().split('\t')
            if lis[0] == '@SQ':
                temp1, ref = lis[1].split(':')
                assert temp1 == 'SN'
                temp2, ref_length = lis[2].split(':')
                assert temp2 == 'LN'
                d_ref_length[ref] = ref_length
                continue
            elif lis[0].startswith('@'):
                #e.g. "@PG"
                continue
            
            flag = int(lis[1])
            qname = lis[0]
            rname = lis[2]
            seq = lis[9]
            pos = lis[3]  # starting aligning position at ref
            if flag & 4:
                # print seq to prefix.nonrna
                print >> fw_nonrna, '>%s\n%s' %(qname, seq)
                continue

            cigar_string = lis[5]    # 
            ### get mismatch
            nm_tag = lis[11]         # NM:i:6, position for bwa-mem
            #nm_tag = lis[?]         # NM:i:6, position for bowtie2 not fixed
            #nm_tag = lis[13]        # NM:i:6, position for bowtie
            for i in lis[11:]:
                if 'NM:i' in i:
                    nm_tag = i
                    break
            nm = int(nm_tag.split(':')[-1])


            read_length = len(seq)


            # get read length in the alignment, M and I in CIGAR
            cigar = parse(cigar_string)
            
            read_align_length = 0
            ref_align_length = 0
            match = 0
            for op, l in cigar:
                if op == 0: # M
                    read_align_length += l
                    ref_align_length += l
                    match += l
                elif op == 1: # I
                    read_align_length += l
                elif op == 2: # D
                    # no base at read
                    # base at ref
                    ref_align_length += l
                    pass
                elif op == 3: # N, many Ds, like intron
                    # same as D
                    ref_align_length += l
                    pass

            coverage = 100.0*read_align_length/read_length
            identity = 100.0*(match - nm)/read_align_length
            read_identity = coverage * identity

            if coverage < args.coverage:
                print >> fw_nonrna, '>%s\n%s' %(qname, seq)
                continue
            if identity < args.align_identity:
                print >> fw_nonrna, '>%s\n%s' %(qname, seq)
                continue
            if read_identity < args.read_identity:
                print >> fw_nonrna, '>%s\n%s' %(qname, seq)
                continue
            if (match - nm) < args.match_number:
                print >> fw_nonrna, '>%s\n%s' %(qname, seq)
                continue

            # qname, n, rname, pos, ref_align_length, read_length, coverage, 
            #  identity, seq
            new_list = [qname, str(n), str(read_length), 
                        rname, str(d_ref_length[rname]), 
                        pos, str(ref_align_length),
                        str(coverage), str(identity), seq]

            print >> fw_temp, '\t'.join(new_list)
            old_read_identity, old_index = d.setdefault(qname, [-1, '-1'])
            if read_identity >= old_read_identity:
                d[qname] = read_identity, str(n)

    d_name = dict((qname, d[qname][-1]) for qname in d)
    with open('%s/%s.rrna' %(args.outdir, args.prefix), 'wb') as fw_rna,\
         open('%s/%s.csv' %(args.outdir, args.prefix), 'wb') as fw_csv:
        ## iter through the temp file to get best hit of multiple hits
        for line in open('%s/%s.temp' %(args.outdir, args.prefix)):
            lis = line.rstrip().split('\t')
            name = lis[0]
            index = lis[1]
            if d_name[name] == index:
                print >> fw_rna, '>%s\n%s' %(lis[0], lis[-1])
                # remove index and seq
                print >> fw_csv, '\t'.join(lis[:1] + lis[2:-1])

if __name__ == '__main__':

    main()
