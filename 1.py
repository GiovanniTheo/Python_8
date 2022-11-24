#!/usr/bin/env python3

import sys

fasta_filename = sys.argv[1]
fasta_filename = 'Python_08.fasta'
fasta_filehandle = open(fasta_filename, 'r')
seq_comp = dict()
seq_id = None

for line in fasta_filehandle:
    line = line.rstrip()
    if line.startswith('>'):
        seq_def = line.lstrip('>').split(' ', maxsplit=1)
        seq_id = seq_def[0]
        seq_comp[seq_id] = {
            'A': 0, 'T': 0, 'G': 0, 'C': 0}    
    else:
        for nucleotide in line:
            seq_comp[seq_id][nucleotide] += 1
            print(seq_comp)

print("seqName\tA_count\tT_count\tG_count\tC_count")

for seq_id in seq_comp:
    print('\t'.join((seq_id,
        str(seq_comp[seq_id]['A']),
        str(seq_comp[seq_id]['C']),
        str(seq_comp[seq_id]['G']),
        str(seq_comp[seq_id]['T'])
    )))
