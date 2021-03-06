#function to read fasta and create a list of sequences
print("usage:\n")
print("python ECquery.py -e eclist.txt -q EL2.ko.txt -f EL2.fasta\npython ECquery.py -p ath00900 -q EL2.ko.txt -f EL2.fasta\n")
print("This script queries EC numbers and returns a fasta of hits.")
print("The input is a list of EC numbers, the query.co file generated by KAAS, and ENZYMES.txt provided with the script.")
print("Please, enter the EC list (-e) or pathway (-p) and the fasta file to fish for sequences (-f)")
print("Query file name (-q) is optional, default is query.ko.txt")


import argparse
import os
from Bio import SeqIO
import pandas as pd
import pickle

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-e', '--EClist', help='FASTA eclist', default='eclist.txt')
parser.add_argument('-q', '--queryKO', help='query.co filename', default='query.ko.txt')
parser.add_argument('-p', '--pathway', help='pathway', default='none')
parser.add_argument('-f', '--fasta', help='Accession rename key file', required=True)

args = parser.parse_args()

eclist = args.EClist
fasta = args.fasta
queryKO = args.queryKO
pathway = args.pathway

print("EC list file: %s; sequences will be extracted from: %s" %(eclist, fasta))

infile = open('ENZYMES.txt')
lines = infile.read()
infile.close()

ENZYMES = lines.split('\n')

ec_d = {}
ECdef_d = {}

for line in ENZYMES:
	ko = line.split()[0][3:]
	definition = ' '.join(line.split()[1:])
	if '[EC:' in line:
		EC = line.split('[EC:')[1]
		definition = definition.split('[EC:')[0]
		if ' ' in EC:
			moreEC = EC.split()
			for item in moreEC:
				if item in ec_d:
					ec_d[item].append(ko)
					ECdef_d[item].append(definition)
				else:
					ec_d[item]=[ko]
					ECdef_d[item]=[definition]
	else:
		EC = line.split()[1]
		EC = EC.split(';')[0]
	if ' ' in EC:
		pass
	else: 
		if EC in ec_d:
			ec_d[EC].append(ko)
			ECdef_d[EC].append(definition)
		else:
			ec_d[EC]=[ko]
			ECdef_d[EC]=[definition]
"""
pickle.dump(ec_d, open('ec_d.txt', 'wb'))
#d = pickle.load(open('ec_d.txt', 'rb'))

pickle.dump(ECdef_d, open('ECdef_d.txt', 'wb'))

"""
infile = open(queryKO)
inquery = infile.read()
infile.close()

kos = inquery.split('\n')
kos_d = {}
for item in kos:
	if '\t' in item:
		seqname = item.split('\t')[0]
		ko = item.split('\t')[1]
		if ko in kos_d:
			kos_d[ko].append(seqname)
		else:
			kos_d[ko]=[seqname]
	else:
		pass

seq_d = {}
infasta = SeqIO.parse(fasta, 'fasta')
for seq in infasta:
	seq_d[seq.name] = seq.seq

infile = open(eclist)
searchedEC = infile.read()
infile.close()

ECs = searchedEC.split()

foundkos = []
foundnames = []
out = open('outfile.fasta', 'w')

for item in ECs:
	enzymedef = ECdef_d[item][0]
	for ortholog in ec_d[item]:
		if ortholog in kos_d:
			foundnames = kos_d[ortholog]
			for name in foundnames:
				sequence = seq_d[name]
				out.write(">%s %s_%s\n%s\n" % (name, item, enzymedef, sequence))		
		else:
			pass

out.close()
