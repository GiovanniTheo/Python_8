#!/usr/bin/env python3

import sys
#abertura e importação do arquivo fasta
fasta_filename = sys.argv[1]
fasta_filehandle = open(fasta_filename, 'r') 

seq_id = None
seq_string = ''
#criação das variáveis se_id e seq_string
for line in fasta_filehandle:
    line = line.rstrip()

    if line.startswith('>'):
        if len(seq_string) > 0:
            print(seq_id + "-frame" + str(1) + "-codons")
            sep = ''  
            for offset in range(0, len(seq_string), 3):
                print("{0}{1}".format(sep, seq_string[offset:offset+3]), end='')
                sep=' '
            print()
#imprime os sq e codons criados
        seq_def = line.lstrip('>').split(' ', maxsplit=1)
        seq_id = seq_def[0]
        seq_string = ''
        #ficca vazio de novo para uma próxima sequência
    else:
        seq_string += line.upper()

if len(seq_string) > 0:
    print(seq_id + "-frame" + str(1) + "-codons")
    sep = ''  
    for offset in range(0, len(seq_string), 3):
        print("{0}{1}".format(sep, seq_string[offset:offset+3]), end='')
        sep=' '
                
    print()
