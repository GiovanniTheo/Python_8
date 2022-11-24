#!/usr/bin/env python3

import sys
#escrevendo para a a criação do arquivo e aplicando testes de erro caso haja entrada errada

try:
    fasta_filename = sys.argv[1]
    if not (fasta_filename.endswith('.fa') or
            fasta_filename.endswith('.nt') or
            fasta_filename.endswith('.fasta')):
        raise Exception("Filename does not end in either .nt, .fa, nor .fasta")

    fasta_filehandle = open(fasta_filename, 'r') 

except Exception as error:  
    print("Error: " + str(error), file=sys.stderr)
    sys.exit(1)

seq_id = None
seq_string = ''
#criando variaveis com os nucleotideos e usando o strip para dividir a sequencia
valid_nucleotides = set('ATCGN')

for line in fasta_filehandle:
    line = line.rstrip()

    if line.startswith('>'):
        if len(seq_string) > 0:

            for frame in range(3):
                print(seq_id + "-frame" + str(frame+1) + "-codons")
                sep = ''  
                for offset in range(0, len(seq_string), 3):
                    if offset+frame+3 <= len(seq_string):
                        print("{0}{1}".format(sep, seq_string[offset+frame:offset+frame+3]), end='')
                        sep=' '
                print()
                
        seq_def = line.lstrip('>').split(' ', maxsplit=1)
        seq_id = seq_def[0]
        seq_string = ''
        
    else:
        line = line.upper()
        
        for nucleotide in line:
            if nucleotide not in valid_nucleotides:
                raise Exception("Invalid nucleotide character in {0}: {1}".format(
                    seq_id, nucleotide
                ))

        seq_string += line
        
if len(seq_string) > 0:
    for frame in range(3):
        print(seq_id + "-frame" + str(frame+1) + "-codons")
        sep = ''  
        for offset in range(0, len(seq_string), 3):
            if offset+frame+3 <= len(seq_string):
                print("{0}{1}".format(sep, seq_string[offset+frame:offset+frame+3]), end='')
                sep=' '
            
        print()
