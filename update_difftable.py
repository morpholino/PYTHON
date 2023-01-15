import os,re
import pandas as pd
import numpy as np
from Bio import SeqIO,Entrez
#http://etetoolkit.org/docs/2.3/tutorial/tutorial_ncbitaxonomy.html
Entrez.email = ''
Entrez.api_key = ''

#ncbi.update_taxonomy_database()
#known issues: 
#does not store already found accessions properly - always starts from the beginning
#for duplicate LOCs, it will attempt to query the accessions on NCBI anyway

def read_blast_results(filetype):
	seq_d = {}
	#run NT_extract.py to obtain these fastas!
	fastas = [x for x in os.listdir(".") if x.endswith("fasta")]
	for file in fastas:
		for seq in SeqIO.parse(file, "fasta"):
			seqname = parse_sseqid(seq.name)
			seq_d[seqname] = seq.seq

	blast_d = {}
	files = [x for x in os.listdir(".") if x.endswith(filetype)]
	for file in files:
		with open(file) as f:
			for l in f:
				l = l.strip().split("\t")
				query = l[0]
				target = parse_sseqid(l[1])
				pident = float(l[2])
				if pident < 95:
					continue
				try:
					blast_d[query] = {}
					blast_d[query]["target"] = target
					blast_d[query]["sequence"] = seq_d[target]
				except KeyError:
					#fasta sequence missing
					pass

	return blast_d

def get_nucl(accession):
	nucl = Entrez.efetch(db='nucleotide', id=accession, rettype='fasta', retmode='text')
	nucl_record = SeqIO.read(nucl, 'fasta')
	acc = nucl_record.name
	desc = nucl_record.description.replace(acc + " ", "")
	sequence = nucl_record.seq
	print("SEQ retrieved:", acc)

	return desc, sequence

def get_aa(accession):
	handle = Entrez.efetch(db='protein', id=accession, rettype='fasta', retmode='text')
	record = SeqIO.read(handle, 'fasta') #only one record is expected
	#nucl_sequence = record['GBSeq_sequence']
	acc = record.name
	desc = record.description.replace(acc + " ", "")
	sequence = record.seq
	print("SEQ retrieved:", acc)

	return desc, sequence

def get_aa_from_nucl(accession):
	#explanation of db structure:
	"""nucl = Entrez.einfo(db='nucleotide')
	record = Entrez.read(nucl)
	record.keys()
	record["DbInfo"]"""

	nucl = Entrez.efetch(db='nucleotide', id=accession, rettype='gb', retmode='xml')
	record = Entrez.read(nucl)[0] #only one record is expected
	#nucl_sequence = record['GBSeq_sequence']
	features = record['GBSeq_feature-table']
	features = [x for x in features if x['GBFeature_key'] == 'CDS'][0] #only one CDS will be extracted
	#orgn = nucl_record.annotations['organism']
	for item in features['GBFeature_quals']:
		if item['GBQualifier_name'] == 'protein_id':
			prot_id = item['GBQualifier_value']
		if item['GBQualifier_name'] == 'translation':
			prot_sequence = item['GBQualifier_value']
			break
	print("SEQ retrieved:", accession)
	sequence = prot_id, prot_sequence #if CDS should be translated

	return sequence

def get_data_from_gene(accession):
	#to retrieve the protein id, the easiest is to parse the gene table
	handle = Entrez.efetch(db='gene', id=accession, rettype='gene_table', retmode='text')
	refseqpattern = r"[A-Z]+_\d+\.\d"
	gene_table = handle.read()
	gene_table = gene_table.split("\n")
	desc = gene_table[0].replace("{} ".format(accession), "")
	desc = desc.replace("[Cyprinus carpio]", "")
	transcript_id = ""
	prot_id = ""
	for line in gene_table:
		if any(x in line[:10] for x in ["RNA"]):
			#>>>>>>>
			#using re, only first isoform is recovered! i.e. LOC109098045
			#caution! loci that are not assigned mRNA annotation! LOC109047210
			try:
				transcript_id = re.findall(refseqpattern, line)[0]
				print("ACC retrieved:", accession, transcript_id)
			except IndexError:
				print("Failed to find accession for RNA:", line)
				continue
			if transcript_id.startswith("XR_"):
				break
		elif line.startswith("protein"):
			#using re:
			prot_id = re.findall(refseqpattern, line)[0] #
			#prot_id = [x for x in line.split() if x.startswith("XP")][0].replace(",", "")
			#prot_id = line.split()[1].replace(",", "")
			print("ACC retrieved:", accession, prot_id)
			if not prot_id.startswith("XP"):
				print(line)
			break

	if prot_id != "":
		#handle = Entrez.efetch(db='protein', id=prot_id, rettype='fasta', retmode='text')
		#record = SeqIO.read(handle, 'fasta')
		#sequence = record.seq
		sequence = ""
		#print("SEQ retrieved:", prot_id)
	else:
		sequence = ""

	if transcript_id == "":
		print("ACC not retrieved:", accession)

	#for xml - more data available:
	#explanation of db structure:
	"""nucl = Entrez.einfo(db='gene')
	record = Entrez.read(nucl)
	record.keys()
	record["DbInfo"]"""
	#handle = Entrez.efetch(db='gene', id=accession, rettype='gb', retmode='xml')
	#record = Entrez.read(handle)[0]
	#print(record.keys())
	#sequence accessions are in ['Entrezgene_comments']

	return transcript_id, prot_id, desc, sequence

def parse_sseqid(sseqid):
	if sseqid.startswith("gi|") or sseqid.startswith("gb|"):
		accession = sseqid.split("|")[1]
	else:
		accession = sseqid
	return accession

##############
###  MAIN  ###
##############
filetype = "met-gtf.tsv"
files = [x for x in os.listdir(".") if x.endswith(filetype)]
rint(files)

blast_d = read_blast_results("blastn")
#print(blast_d)

for file in files:
	print("\n\n======\nAnalyzing", file)
	accessionset = set()
	protset = set()
	processedfile = file.replace(filetype, "oook")
	if os.path.exists(processedfile):
		print("record file exists")
		with open(processedfile) as f:
			for l in f:
				l = l.strip().split()
				accessionset.add(l[3])
				if l[3].startswith("XP"):
					protset.add("{}_{}-{}".format(l[0], l[1], l[2]))
			print(len(accessionset), "items finished")
	dataset = file.split(".")[0]
	processed = open(processedfile, "a")
	with open(file) as f:
		records = pd.DataFrame(columns=["Chr","Chr.start","Chr.end",
										"TXX.Mean.Diff","pDMR","qDMR.BH",
										"seqnames","start","end","width","strand","type",
										"gene_id","transcript_id","gene_biotype",
										"product","transcript_biotype",
										"protein_id","description"])
		for c,l in enumerate(f):
			c += 1
			#i need the following data: chromosomal coordinates l[0:3], LOC* locus ID l[12], XM* mRNA ID l[15], 
			#"gene_type" l[14], XP* protein ID l[17], annotation l[16]
			line = l.strip().split("\t")
			if line[-2] == "protein_id":
				continue

			#query = "{}_{}-{}".format(line[0], line[1], line[2]) #chromosome
			records.loc[c, "Chr"] = line[0]
			records.loc[c, "Chr.start"] = line[1]
			records.loc[c, "Chr.end"] = line[2]
			records.loc[c, "TXX.Mean.Diff"] = line[3]
			records.loc[c, "pDMR"] = line[4]
			records.loc[c, "qDMR.BH"] = line[5]
			records.loc[c, "seqnames"] = line[6]
			records.loc[c, "start"] = line[7]
			records.loc[c, "end"] = line[8]
			records.loc[c, "width"] = line[9]
			records.loc[c, "strand"] = line[10]
			records.loc[c, "type"] = line[11]
			records.loc[c, "product"] = line[15]
			records.loc[c, "transcript_biotype"] = line[16]


			aa_accession = line[17]
			if aa_accession.startswith("XP"):
				records.loc[c, "protein_id"] = aa_accession
				description,_ = get_aa(aa_accession)
				records.loc[c, "description"] = description
				#accessionset.add(aa_accession) #only after the fetch, accession is added
				#protset.add(aa_accession)
				processed.write("{}_P: {}\t{}\t{}\t{}\n".format(c, line[0], line[1], line[2], aa_accession))
				processed.write("{}_D: {}\t{}\t{}\t{}\n".format(c, line[0], line[1], line[2], description))
			elif aa_accession != "NA":
				records.loc[c, "protein_id"] = aa_accession

			transcript_id = line[13]					
			if transcript_id.startswith("X"):
				records.loc[c, "transcript_id"] = transcript_id
				description,_ = get_nucl(transcript_id)
				if records.loc[c, "description"] in ["", "NA", np.nan]:
					records.loc[c, "description"] = description
				processed.write("{}_T: {}\t{}\t{}\t{}\n".format(c, line[0], line[1], line[2], transcript_id))
				processed.write("{}_D: {}\t{}\t{}\t{}\n".format(c, line[0], line[1], line[2], description))
			#elif query in blast_d:
			#	sequence = blast_d[query]["sequence"]
			#	accession = blast_d[query]["target"]

			gene_id = line[12]
			if gene_id.startswith("LOC"):
				records.loc[c, "gene_id"] = gene_id
				#try finding protein ID from locus info
				if aa_accession in ["", "NA", np.nan]:
					transcript_id, aa_accession, description, _ = get_data_from_gene(gene_id)
					#retry aa_accession retrieval
					if aa_accession == "":
						print("problem retrieving protein ID", gene_id)
					else:
						records.loc[c, "protein_id"] = aa_accession
						processed.write("{}_L: {}\t{}\t{}\t{}\t{}\n".format(c, line[0], line[1], line[2], gene_id, aa_accession))
					if records.loc[c, "transcript_id"] in ["", "NA", np.nan]:
						records.loc[c, "transcript_id"] = transcript_id
						processed.write("{}_L: {}\t{}\t{}\t{}\t{}\n".format(c, line[0], line[1], line[2], gene_id, transcript_id))
					elif transcript_id == "":
						print("problem retrieving transcript ID", gene_id)
					else:
						#transcript ID already stored in dataframe
						pass
					if records.loc[c, "description"] in ["", "NA", np.nan]:
						records.loc[c, "description"] = description
			elif gene_id != "NA":
				records.loc[c, "gene_id"] = gene_id
			

			gene_biotype = line[14]
			if gene_biotype != "NA":
				records.loc[c, "gene_biotype"] = gene_biotype
			elif records.loc[c, "protein_id"] not in ["", "NA", np.nan]:
				records.loc[c, "gene_biotype"] = "protein_coding"
			elif records.loc[c, "transcript_id"] == "XR":
				records.loc[c, "gene_biotype"] = "lncRNA"

			description = line[18]
			if description != "NA":
				#override if annotation has been provided
				records.loc[c, "description"] = description

	processed.close()
	records.to_csv(file.replace("tsv", "new.tsv"), sep="\t")
			
print("Finished")

