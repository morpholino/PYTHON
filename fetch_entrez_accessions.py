#!/usr/bin/python3
import os
import re
from Bio import Entrez,SeqIO
from ete3 import NCBITaxa
#http://etetoolkit.org/docs/2.3/tutorial/tutorial_ncbitaxonomy.html
ncbi = NCBITaxa()

Entrez.email = 'zoltan.fussy@gmail.com'
Entrez.api_key = "ed51bca6c73792ecc692af11ff762b72a008"

if os.path.isdir("/Users/morpholino/OwnCloud/"):
	home = "/Users/morpholino/OwnCloud/"
elif os.path.isdir("/Volumes/zoliq data/OwnCloud/"):
	home = "/Volumes/zoliq data/OwnCloud/"
else:
	print("Please set a homedir")

#wd = "AndyLab/Lysine pathway/v2"
wd = "AndyLab/Phaeocystis/annotation/polyamines"
#wd = "DolezalLab/SecY/"
os.chdir(home + wd)

#will need these regexes
taxonpattern = r'\[(.+)\]'
genbankpattern = r'[A-Z]{1,3}(_)?\d+(\.\d)?'

#collect files
fastas = [x for x in os.listdir(".") if x.endswith("_out.fasta")] #NOTE: default output suffix of the diamondparse script
fastas.sort()

#main
for file in fastas:
	print("=============\n\n\nNow parsing file", file)
	fasta = SeqIO.parse(file, "fasta")
	with open(file.replace("_out.fasta", "_v1.fasta"), 'w') as out:
		for f in fasta:
			fname = f.name.strip()
			prot_idregex = re.search(genbankpattern, fname)
			#print(fname)
			if prot_idregex:
				prot_id = prot_idregex.group(0)
				#print("GenBank ID found, yay!", prot_id)

			#https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly
				if fname.startswith(prot_id):
					print("Fetching data", prot_id)
					try:
						prot = Entrez.efetch(db='protein', id=prot_id, rettype='fasta', retmode='text')
						prot_record = SeqIO.read(prot, 'fasta')
						desc = prot_record.description
						taxon = re.search(taxonpattern, desc).group(1)
						if taxon:
							taxon = "_".join(taxon.split()[:2])
							desc = f'{taxon}_{desc.split(" [")[0]}'
						seq = prot_record.seq
						print("success", desc, seq[:20], "...")
					except:
						#bad request
						desc = f.name
						seq = f.seq
				else:
					desc = f.description
					seq = f.seq
			else:
				desc = f.description
				seq = f.seq
			out.write(f">{desc}\n{seq}\n")
			#out.write('{}\t{}__{}\n'.format(prot_id, orgn, tax))


quit("\n\nFinished processing, yay!")
#this is rich format:
for prot_id in ids:
	#print(prot_id)
	prot = Entrez.efetch(db='protein', id=prot_id, rettype='gb', retmode='text')
	prot_record = SeqIO.read(prot, 'genbank')
	print(prot_record)
	#tax = str(prot_record.annotations['taxonomy'][::-1]).replace('\'', '')
	orgn = prot_record.annotations['organism']
	name2taxid = ncbi.get_name_translator(orgn)
	print(prot_id, orgn, name2taxid)
	#out.write('{}\t{}__{}\n'.format(prot_id, orgn, tax))