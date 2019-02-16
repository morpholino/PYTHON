import os
import csv
from Bio import Entrez,SeqIO
Entrez.email = "zoltan.fussy@google.com" 

def getseq(accession):
	try:
		handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="XML")
		record = Entrez.read(handle)
		try:
			sequence = record[0]['TSeq_sequence']
		except KeyError:
			sequence = "-"
	except:
		print("error retrieving", accession)
		sequence = "-"
	return sequence

############
### MAIN ###	
############

#pick what sort of files you need to parse:
filetype =  ["fasta", "renamefile", "tsv", "errorfile-missing"][2]

badchars = ("|@+,:;()'")
allowed = ("fasta", "fas", "fst", "phy", "phylip")

home = "/Users/zoliq/ownCloud/"
#home = "/Volumes/zoliq data/OwnCloud/"
#wd = home + "Jankoviny/Tick_transcriptome/lectintree-filt"
wd = home + "genomes/chromera/plastid proteome/"
os.chdir(wd)

if filetype == "renamefile":
	files = [x for x in os.listdir(".") if x.startswith("rename")]
	for file in files:
		newfile = "new_" + file
		print("Processing file " + file)
		with open(file) as f, open(newfile, "w") as out:
			for l in f:
				if l.startswith("gi_"):
					gid = l.split("_")[1]
					handle = Entrez.efetch(db="protein", id=gid, rettype="fasta", retmode="XML")
					record = Entrez.read(handle)
					try:
						accession = record[0]['TSeq_accver']
					except KeyError:
						accession = l.split("\t")[0].replace("gi_{}_".format(gid), "")
					#print(record)
					annot = "".join([x for x in record[0]['TSeq_defline'] if x not in badchars])
					#print(annot)
					out.write("{}\t{}_{}\n".format(l.split("\t")[0], accession, annot))
					#record['TSeq_taxid'] and record['TSeq_orgname'] also possible
				else:
					out.write(l)

elif filetype == "fasta":
	files = [x for x in os.listdir(".") if x.endswith("fasta")]
	#test purposes only: files = ["STT3_pfam02516.fasta"]
	for file in files:
		newfile = "new_" + file
		print("Processing file " + file)
		f = SeqIO.parse(file, "fasta")
		with open(newfile, "w") as out:
			for l in f:
				if l.name.startswith("gi|"):
					gid = l.name.split("|")[1]
					handle = Entrez.efetch(db="protein", id=gid, rettype="fasta", retmode="XML")
					record = Entrez.read(handle)
					try:
						accession = record[0]['TSeq_accver']
					except KeyError:
						accession = l.name.replace("gi|{}|".format(gid), "")
					#print(record)
					annot = "".join([x for x in record[0]['TSeq_defline'] if x not in badchars])
					#print(annot)
					out.write(">{} {}\n{}\n".format(accession, annot, l.seq))
					#record['TSeq_taxid'] and record['TSeq_orgname'] also possible
				else:
					description = "".join(x for x in l.description if x not in badchars)
					out.write(">{}\n{}\n".format(description, l.seq))

#this part is to retrieve accessions written in a table
elif filetype == "tsv":
	files = [x for x in os.listdir(".") if x.endswith("tsv")]
	#test purposes only: files = ["STT3_pfam02516.fasta"]
	for file in files:
		missing = file.replace(".tsv", "") + "-missing.txt"
		newfile = file.replace(".tsv", "") + "-get.fasta"
		print("Processing file " + file)
		f = csv.reader(open(file), delimiter="\t", skipinitialspace=False)
		with open(newfile, "w") as out, open(missing, "w") as error:
			for c,l in enumerate(f):
				print(c, l)
				if c == 0:
					species = l.copy()
					continue
				else:
					name = l[0]
					accessions = l[1:]
				for i,a in enumerate(accessions):
					tag = "{}_{}".format(name, species[i+1])
					if "@" in a:
						borders = a.split("@")[0].split("-")
						queries = [a.split("@")[1]]
						print("sequence splitting not finished")
					elif "&" in a:
						queries = a.split("&")
						for q in queries:
							seq = getseq(q)
							if seq != "-":
								out.write(">{}_{}\n{}\n".format(q, tag, seq))
							else:
								error.write("{}\t{}\n".format(q, tag))
					elif a in ["", "-"]:
						pass
					else:
						seq = getseq(a)
						if seq != "-":
							out.write(">{}_{}\n{}\n".format(a, tag, seq))
						else:
							error.write("{}\t{}\n".format(a, tag))
						print("getseq function on single accession")

#this part is to search for missed sequences:
elif filetype == "errorfile-missing":
	files = [x for x in os.listdir(".") if x.endswith("-missing.txt")]
	for file in files:
		missingdict = {}
		missing = file.replace("missing", "missing2")
		newfile = file.replace(".txt", "") + "-err_get.fasta"
		with open(file) as f:
			for l in f:
				data = l.strip().split("\t")
				if len(data) == 2:
					missingdict[data[0]] = data[1]
		with open(missing, "w") as error, open(newfile, "w") as result:
			infasta1 = SeqIO.parse(home + "genomes/sources/euks/PlasmoDB-40_Pfalciparum3D7_AnnotatedProteins.fasta", "fasta")
			for seq in infasta1:
				name = seq.name.split(".")[0]
				if name in missingdict:
					tag = missingdict[name]
					result.write(">{}_{} {}\n{}\n".format(name, tag, seq.description, seq.seq))
					del missingdict[name]

			infasta2 = SeqIO.parse(home + "genomes/sources/euks/ToxoDB-40_TgondiiME49_AnnotatedProteins.fasta", "fasta")
			for seq in infasta2:
				name = seq.name.split("-")[0]
				if name in missingdict:
					tag = missingdict[name]
					result.write(">{}_{} {}\n{}\n".format(name, tag, seq.description, seq.seq))
					del missingdict[name]
					altname = name.replace("TGME49_2", "TGME49_0").replace("TGME49_3", "TGME49_1")
					try:
						del missingdict[altname]
					except KeyError:
						print(altname, "not in dict")
			
			infasta3 = SeqIO.parse(home + "genomes/sources/euks/cyanidioschyzon_aa.fasta", "fasta")
			for seq in infasta3:
				name = seq.name.split(".")[0]
				if name in missingdict:
					tag = missingdict[name]
					result.write(">{}_{} {}\n{}\n".format(name, tag, seq.description, seq.seq))
					del missingdict[name]
			
			for left in missingdict:
				error.write("{}\t{}\n".format(left, missingdict[left]))
"""
Swissprot works in a similar way:
>>> from Bio import ExPASy,SwissProt

>>> handle = ExPASy.get_sprot_raw(hitid)
>>> record = SwissProt.read(handle)
>>> dir(record)
['__doc__', '__init__', '__module__', 'accessions', 'annotation_update',
'comments', 'created', 'cross_references', 'data_class', 'description',
'entry_name', 'features', 'gene_name', 'host_organism', 'keywords',
'molecule_type', 'organelle', 'organism', 'organism_classification',
'references', 'seqinfo', 'sequence', 'sequence_length',
'sequence_update', 'taxonomy_id']
"""