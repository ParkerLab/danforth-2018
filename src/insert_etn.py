#!/usr/bin/env python

#
# Inserts the retrotransposon sequence at its position in mm9 chromosome 2.
#

LINE_LENGTH = 50

etn = ''
with open('JX863104.fasta', 'r') as etn_file:
    next(etn_file)  # skip header line
    etn = ''.join((line.strip() for line in etn_file))

print('Length of etn: ', len(etn))

chr2_id = ''
chr2 = ''
with open('mm9.chr2.fasta') as chr2_file:
    chr2_id = next(chr2_file)
    chr2 = ''.join(line.strip() for line in chr2_file)

print('Original chr2 length: ', len(chr2))

chr2etn = chr2[:19355026] + etn[6:] + chr2[19355026:]
chr2etn_len = len(chr2etn)
print('Modified chr2 length: ', chr2etn_len)

with open('mm9.chr2.etn.fasta', 'w') as chr2etn_file:
    chr2etn_file.write(chr2_id)
    i = 0
    while i < chr2etn_len:
        chr2etn_file.write(chr2etn[i:i + LINE_LENGTH])
        chr2etn_file.write('\n')
        i += LINE_LENGTH
    chr2etn_file.close()
