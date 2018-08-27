#function to read fasta and create a list of sequences
print("usage:\n")
print("python ECquery.py -e eclist.txt -q query.ko.txt -f query.fasta\nor\npython ECquery.py -p ko00900 -q query.ko.txt -f query.fasta\n")
print("This script queries EC numbers (accepts also KEGG ortholog numbers) and returns an annotated fasta of hits.")
print("The input is a list of EC numbers, the query.co file generated by KAAS, and ENZYMES.txt provided/prepared with the script.")
print("Please, enter the EC list (-e) or pathway (-p) and the fasta file to fish for annotations (-f)")
print("Query file name (-q) is optional, default is query.ko.txt")


import argparse
import os
import urllib
from Bio import SeqIO
#import pandas as pd
#import pickle
from bioservices.kegg import KEGG

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-e', '--EClist', help='EC list to download', default='eclist.txt')
parser.add_argument('-q', '--queryKO', help='query.ko filename', default='query.ko.txt')
parser.add_argument('-p', '--pathway', help='pathway to download (optional)', default='none')
parser.add_argument('-f', '--fasta', help='queried FASTA file', required=True)

args = parser.parse_args()

eclist = args.EClist
fasta = args.fasta
queryKO = args.queryKO
pathway = args.pathway

print("EC list file: %s; sequences will be extracted from: %s" %(eclist, fasta))

# UNCOMMENT TO DOWNLOAD THE CURRENT LIST OF KEGG ENZYMES FROM THE DATABASE
# MIGHT HELP WHEN GETTING KEY ERRORS OF THE EC FILE
"""
url = 'http://rest.kegg.jp/list/ko'
fileName = 'ENZYMES.txt'
urllib.request.urlretrieve(url, fileName)
print("Current list of enzymes downloaded.")
"""
#first, dictionaries of ECs and their ortholog groups is made
#ortholog groups are called from ec_d, definitions from ECdef_d
#retrieve up to date list visiting: http://rest.kegg.jp/list/ko
infile = open('ENZYMES.txt')
lines = infile.read()
infile.close()


ENZYMES = lines.split('\n')[:-1]


KO_dic = {}
EC2KO_dic = {}

for line in ENZYMES:
	ko = line.split()[0][3:]
	#print(ko)
	definition = ' '.join(line.split()[1:])
	if '[EC:' in definition:
		EC = line.split('[EC:')[1]
		EC = EC.replace("]","")
		definition = definition.split('[EC:')[0]
		if ' ' in EC:
			EC = EC.split()
		else:
			EC = [EC]
	
		KO_dic[ko] = [definition, EC]
		for item in EC:
			if item in EC2KO_dic:
				EC2KO_dic[item].append(ko) #here append may be used, bc only one item is added at a time
			else:
				EC2KO_dic[item] = [ko]
	else:
		KO_dic[ko] = [definition]


#here a list of KOs is made to be searched for - if pathway is defined
kostoget = {}

if args.pathway != 'none':
	output = open('eclist.txt', 'w')

	kegg = KEGG()
	pathway = kegg.get(args.pathway)
	dict_data = kegg.parse(pathway)
	#print(dict_data)
	try:
		for key in dict_data['ORTHOLOGY'].keys():
			print("adding ortholog {} to eclist".format(key))
			value = dict_data['ORTHOLOGY'][key]
			print(value)
			if "," in key:
				for item in key.split(","):
					kostoget[item] = value
					output.write(item + '\n')
			else:
				kostoget[key] = value 
				output.write(key + '\n')
		print("eclist.txt written")
	except KeyError as modularpathway:
		modules = []
		for key in dict_data['MODULE'].keys():
			modules.append(key)
		print("Modular pathway selected, do search for: " + ", ".join(modules))
		for module in modules:
			pathway = kegg.get(module)
			dict_data = kegg.parse(pathway)
			for key in dict_data['ORTHOLOGY'].keys():
				print("adding ortholog {} to eclist".format(key))
				value = dict_data['ORTHOLOGY'][key]
				print(value)
				if "," in key:
					for item in key.split(","):
						kostoget[item] = value
						output.write(item + '\n')
				else:
					kostoget[key] = value 
					output.write(key + '\n')

		print(kostoget)
	output.close()

#here query.ko.txt is opened and for each orthology group respective contigs are saved to kos_d
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
print("queryKO file analyzed")

#fasta file is written and contigs are saved to seq_d
seq_d = {}
infasta = SeqIO.parse(fasta, 'fasta')
for seq in infasta:
	seq_d[seq.name] = seq.seq


infile = open(eclist)
ECs = infile.read().split()
infile.close()
kostoprint = set() #this was changed from list to set() so that KOs do not repeat
for item in ECs:
	if item.startswith("K"):
		kostoprint.append(item)
	else:
		try:
			kostoprint |= set(EC2KO_dic[item])
			#here += must be used for lists, otherwise lists of lists will be created
			#in case of sets, use union: kostoprint |= set(EC2KO_dic[item])
		except KeyError:
			print(item + "not in list of ECs")
#print(kostoprint)

out = open('outfile.fasta', 'w')

foundnames = []
for ortholog in kostoprint:
	if ortholog in kos_d.keys():
		foundnames = kos_d[ortholog]
		enzymedef = ortholog

		#to print enzyme definitions to fasta desription
		if len(KO_dic[ortholog]) == 1:	
			enzymedef = KO_dic[ortholog][0]
		elif len(KO_dic[ortholog]) == 2:
			twopartdef = [KO_dic[ortholog][0], ", ".join(KO_dic[ortholog][1])]
			enzymedef = " ".join(twopartdef)

		for name in foundnames:
			sequence = seq_d[name]
			out.write(">{}_{} {}\n{}\n".format(name, ortholog, enzymedef, sequence))
	else:
		pass
		#print("Ortholog {} not in the given dataset.".format(ortholog))


for ortholog in kostoget: 
	if ortholog in kos_d:
		foundnames = kos_d[ortholog]
		for name in foundnames:
			sequence = seq_d[name]
			out.write(">{} {}\n{}\n".format(name, KO_dic[ortholog], sequence))		
	else:
		pass		
		#print("Ortholog {} not in the given dataset.".format(ortholog))

out.close()

print("now dereplicating...")

output1 = open('outfile_dedupl_desc.fasta', 'w')
output2 = open('outfile_dupl_names.fasta', 'w')

multiplications = {}
seq_dict = {}
for sequence in SeqIO.parse('outfile.fasta', 'fasta'):
	if sequence.seq in multiplications:
		multiplications[sequence.seq].append(sequence.name)
	else:
		multiplications[sequence.seq] = [sequence.name]
	if sequence.seq not in seq_dict:
		#rename full header only with name (acc, till the first space)
		# seq_dict[sequence.seq] = sequence.name 
		#keep full header			
		seq_dict[sequence.seq] = sequence.description

for key, value in seq_dict.items():
	output1.write('>{}\n{}\n'.format(value, key))

for key, value in multiplications.items():
	if len(value) > 1:
		output2.write('{}\n'.format(str(value)))

output1.close()
output2.close()

print("dereplicated hits written to outfile_dedupl_desc.fasta")
