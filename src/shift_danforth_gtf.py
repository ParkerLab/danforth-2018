#!/usr/bin/env python3

#
# shift_danforth: adjusts the mm9 GTF annotation file for the Danforth early transposon retroviral insertion
#

#
# zcat /lab/data/reference/mouse/mm9/annot/gencode.vM1.annotation.gtf.gz | ./shift_danforth -v 2> shifts.txt| gzip > annotation.gtf.gz
#

import argparse
import errno
import os.path
import signal
import sys

ETN_REFERENCE = 'chr2'
ETN_START = 19355026
ETN_LENGTH = 8534
ETN_END = ETN_START + ETN_LENGTH


def overlaps_etn(start, end):
    return (ETN_START <= start <= ETN_END) or (ETN_START <= end <= ETN_END)


argparser = argparse.ArgumentParser(
    prog=sys.argv[0],
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="""Transforms the given mm9 annotation GTF to shift features after the ETn at {}:{}\n""".format(ETN_REFERENCE, ETN_START)
)

argparser.add_argument('-v', '--verbose', action='store_true', help="""Report all features shifted.""")
argparser.add_argument('mm9', nargs='?', default='-', help="""The GTF file to transform, or '-' to read from standard input.""")
argparser.add_argument('danforth', nargs='?', default='-', help="""The GTF file to create, or '-' to write to standard output.""")

args = argparser.parse_args()

if args.mm9 != '-' and not os.path.isfile(args.mm9):
    print('Cannot open annotation file {}'.format(args.annotation), file=sys.stderr)
    sys.exit(1)

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

danforth = args.danforth == '-'and sys.stdout or open(args.danforth, 'w')
mm9 = args.mm9 == '-' and sys.stdin or open(args.mm9)

for line in mm9:
    if line[0] == '#':
        try:
            danforth.write(line)
        except IOError as e:
            if e.errno == errno.EPIPE and args.mm9 == '-':
                pass
    else:
        seqname, source, feature, start, end, score, strand, frame, attribute = line.split('\t')
        if seqname == ETN_REFERENCE:
            start = int(start)
            end = int(end)
            if ETN_START <= start:
                if args.verbose:
                    newstart = start + ETN_LENGTH
                    newend = end + ETN_LENGTH
                    print('Shifting {}:{}-{}{} to {}:{}-{}'.format(
                        seqname,
                        start,
                        end,
                        overlaps_etn(start, end) and ' (overlaps ETn at {}:{}-{})'.format(ETN_REFERENCE, ETN_START, ETN_END) or '',
                        seqname,
                        newstart,
                        newend
                    ), file=sys.stderr)
                    start = newstart
                    end = newend

        danforth.write('\t'.join(str(s) for s in [seqname, source, feature, start, end, score, strand, frame, attribute]))
