#!/usr/bin/python3
from Bio import SeqIO

#infile = SeqIO.parse('./Trinity_Sterlet_mat,oos_WTA.fa', 'fasta')
#output = open('./Trinity_Sterlet_mat,oos_WTA_6frameAA.fa', 'w')
infile = SeqIO.parse('./ryba_filtfasta.fasta', 'fasta')
output = open('./ryba_filtfasta_6frameAA.fasta', 'w')

gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'X', 'TAG':'X',
    'TGC':'C', 'TGT':'C', 'TGA':'X', 'TGG':'W'}

def translation(sequence):
    cut_seq = []
    for i in range(0,len(sequence)-2,3):
        cut_seq.append(sequence[i:i+3])
    aa = []
    for codon in cut_seq:
        if 'N' in codon:
            aa.append('X')
        else:
            aa.append(gencode[str(codon)]) #with the str() here also runs with python2.7
    return ''.join(aa)

i = 0
for sequence in infile:
    i += 1
    if i%10000 == 0:
        print(i)
    name = sequence.name
    seq = sequence.seq.upper()
    ambiguous = False
    for nucleotide in seq:
        if nucleotide not in 'ATCGN':
            ambiguous = True
    if ambiguous == True:
        with open('error.fasta', 'a') as error:
            error.write('{}\n{}\n'.format(name, seq)) 
    output.write('>{}_1\n{}\n'.format(name, translation(seq)))
    output.write('>{}_2\n{}\n'.format(name, translation(seq[1:])))
    output.write('>{}_3\n{}\n'.format(name, translation(seq[2:])))
    output.write('>{}_4\n{}\n'.format(name, translation(seq.reverse_complement())))
    output.write('>{}_5\n{}\n'.format(name, translation(seq.reverse_complement()[1:])))
    output.write('>{}_6\n{}\n'.format(name, translation(seq.reverse_complement()[2:])))

output.close()
