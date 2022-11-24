#!/usr/bin/env python3

import sys

try:
    fasta_filename = sys.argv[1]

    if not (fasta_filename.endswith('.fa') or
            fasta_filename.endswith('.nt') or
            fasta_filename.endswith('.fasta')):
        raise Exception("Filename does not end in either .nt, .fa, nor .fasta")

except Exception:
    print("Error: No input file or incorrect file suffix", file=sys.stderr)
    sys.exit(1)

try:
    output_prefix = sys.argv[2]
    
except IndexError:
    print("Error: No output file prefix", file=sys.stderr)
    sys.exit(1)

codon_filename = output_prefix + '.codons-6frames.nt'
trans_filename = output_prefix + '.translated.aa'
long_aa_filename = output_prefix + '.translated-longest.aa'
long_nt_filename = output_prefix + '.orf-longest.nt'


try:    
    fasta_file = open(fasta_filename, 'r')
    codon_file = open(codon_filename, 'w')
    trans_file = open(trans_filename, 'w')
    long_aa_file = open(long_aa_filename, 'w')
    long_nt_file = open(long_nt_filename, 'w')
except IOError as error:
    print("Error: " + str(error), file=sys.stderr)
    sys.exit(1)    
    
seq_id = None
seq_string = ''

valid_nucleotides = set('ATCGN')
complement_table = str.maketrans({'A':'T','C':'G','T':'A','G':'C','N':'N'})
translation_table = {
    'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
    'AAT':'N', 'AAC':'N',
    'GAT':'D', 'GAC':'D',
    'TGT':'C', 'TGC':'C',
    'CAA':'Q', 'CAG':'Q',
    'GAA':'E', 'GAG':'E',
    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
    'CAT':'H', 'CAC':'H',
    'ATT':'I', 'ATC':'I', 'ATA':'I',
    'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
    'AAA':'K', 'AAG':'K',
    'ATG':'M',
    'TTT':'F', 'TTC':'F',
    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'TGG':'W',
    'TAT':'Y', 'TAC':'Y',
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
    'TAA':'*', 'TGA':'*', 'TAG':'*'
}
    
codon_length = 3
seq_ids = list()
seq_strings = dict()
for line in fasta_file:
    line = line.rstrip()

    if line.startswith('>'):
        seq_def = line.lstrip('>').split(' ', maxsplit=1)
        seq_id = seq_def[0]
        seq_strings[seq_id] = {'+':'', '-':''}
        seq_ids.append(seq_id)
        
    else:
        line = line.upper()
        for nucleotide in line:
            if nucleotide not in valid_nucleotides:
                raise Exception("Invalid nucleotide character in {0}: {1}".format(
                    seq_id, nucleotide
                ))

        seq_strings[seq_id]['+'] += line

for seq_id in seq_ids:
    seq_strings[seq_id]['-'] = seq_strings[seq_id]['+'].translate(complement_table)

    info = list()
    nt_list = list()
    aa_list = list()
    num_translations = 0
    for strand in ('+', '-'):
        seq_string = seq_strings[seq_id][strand]

        for frame in range(codon_length):
            num_orfs = 0
            codons = []
            aminos = []

            for offset in range(0, len(seq_string), codon_length):
                beg = offset + frame
                end = beg + codon_length

                if end <= len(seq_string):
                    codon = seq_string[beg:end]
                    amino = translation_table[codon]
                                    
                    codons.append(codon)
                    aminos.append(amino)
                    
                    if amino == '*':
                        info.append((frame, num_orfs))
                        nt_list.append(codons)
                        aa_list.append(aminos)
                        num_translations += 1
                        num_orfs += 1
                        codons = []
                        aminos = []
                    
            if len(codons) > 0:
                info.append((frame, num_orfs))
                nt_list.append(codons)
                aa_list.append(aminos)
                num_translations += 1
                num_orfs += 1                
        
    longest_frame = -1
    longest_length = -1
    longest_codons = []
    longest_aminos = []
    for i in range(num_translations):  
        codons = nt_list[i]
        aminos = aa_list[i]
        frame, orf = info[i]

        if aminos[0] == '*':
            continue
        
        frame += (codon_length if info[i] == '-' else 0)  

        try:
            first_M = aminos.index('M')
        except ValueError:
            
            continue

        if aminos[-1] != '*':
            continue

        orf_codons = codons[first_M:]
        orf_aminos = aminos[first_M:]
        
        codon_file.write('>{0}-frame-{1}-orf-{2}-codons\n'.format(seq_id, frame+1, orf+1))
        codon_file.write(' '.join(orf_codons) + '\n')
        trans_file.write('>{0}-frame-{1}-orf-{2}-aminos\n'.format(seq_id, frame+1, orf+1))
        trans_file.write(''.join(orf_aminos) + '\n')

        if len(orf_aminos) > longest_length:
            longest_length = len(orf_aminos)
            longest_codons = orf_codons
            longest_aminos = orf_aminos
            longest_frame  = frame
            longest_orf    = orf
            
    long_aa_file.write('>{0}-frame-{1}-orf-{2}-aminos\n'.format(seq_id, longest_frame+1, longest_orf+1))
    long_aa_file.write(''.join(longest_aminos) + '\n')
    long_nt_file.write('>{0}-frame-{1}-orf-{2}-codons\n'.format(seq_id, longest_frame+1, longest_orf+1))
    long_nt_file.write(' '.join(longest_codons) + '\n')
