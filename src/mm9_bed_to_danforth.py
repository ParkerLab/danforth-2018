#!/usr/bin/env python

import sys

INSERTION_CHROMOSOME = 'chr2'
INSERTION_POSITION = 19355026 - 1
INSERTION_LENGTH = 8534


def usage():
    message = 'Description: shift a bed file from mm9 to danforth coordinates\n'
    message += 'First three columns are presumed to be chrom/start/end\n'
    message += '\nUsage: {} input.bed > input.shifted.bed\n'.format(sys.argv[0])
    sys.stderr.write(message)


if len(sys.argv) == 1:
    usage()
    sys.exit()

bed = sys.argv[1]

with open(bed, 'r') as f:
    for line in f:
        line = line.rstrip()
        line = line.split("\t")

        if line[0] == INSERTION_CHROMOSOME or line[0] == INSERTION_CHROMOSOME.replace('chr', ''):
            if int(line[1]) >= INSERTION_POSITION:
                line[1] = str(int(line[1]) + INSERTION_LENGTH)
                line[2] = str(int(line[2]) + INSERTION_LENGTH)

        line = '\t'.join(line)
        print(line)
