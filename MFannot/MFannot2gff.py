import os,argparse

#to parse arguments listed in the command after the script
parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-i', '--infile', help='MFannot output to process', default="userdata.new.tbl")
parser.add_argument('-o', '--outfile', help='output prefix', default="comboContig1")
parser.add_argument('-d', '--directory', help='working directory', default="here")
args = parser.parse_args()

##################################
#### Create working directory ####
##################################

if os.path.isdir("/Users/morpholino/OwnCloud/"):
	home = "/Users/morpholino/OwnCloud/"
elif os.path.isdir("/Volumes/zoliq data/OwnCloud/"):
	home = "/Volumes/zoliq data/OwnCloud/"
else:
	print("Please set a homedir")
if args.directory == "here":
	print("changing to default directory")
	defdir = "AndyLab/Phaeocystis/reassemblies/organellar"
	wd = home + defdir
	os.chdir(wd)
else:
	os.chdir(args.directory)

codechange = {"Ala": "A", "Arg": "R", "Asn": "N",
"Asp": "D", "Cys": "C", "Glu": "E", 
"Gln": "Q", "Gly": "G", "His": "H",
"Ile": "I", "Leu": "L", "Lys": "K",
"Met": "M", "Phe": "F", "Pro": "P",
"Ser": "S", "Thr": "T", "Trp": "W",
"Tyr": "Y", "Val": "V", "Sup": "W", "Undet": "X"}

if args.infile == "batch":
	infiles = os.listdir('.')
	infiles = [f for f in infiles if f.endswith("tbl")]
	#example of these files in bico/genomeseq/spadesfilt/draft genome/mito_proteome/*
else: 
	infiles = [args.infile]

#####################
####   Modules   ####
#####################

def tblread(file):
	tbl_dic = {}
	with open(file) as f:
		for line in f:
			if line.startswith(">Feature"):
				feature = line.split()[1]
				tbl_dic[feature] = {}
			elif len(line.split("\t")) > 1:
				line = line.split("\t")
				if line[0] != "":
					start = int(line[0].replace(">", "").replace("<", ""))
					stop = int(line[1].replace(">", "").replace("<", ""))
					if start < stop:
						direction = "forward"
					else:
						direction = "reverse"
					modulerange = (min(start, stop), max(start, stop))
				elif line[3] == "gene":
					modulename = line[4].strip()
					if modulename.startswith("trn"):
						modulename = modulename[:4]
					tbl_dic[feature][modulename] = {"Range": modulerange, 
											"Database": "MFannot", 
											"Strand": direction,
											"Description": ""}
				elif line[3] == "product":
					#modulename stays the same
					tbl_dic[feature][modulename].update({"Description": line[4].strip()})
	return tbl_dic

##################
####   Main   ####
##################

#infiles = [x for x in infiles if x.endswith("tbl")]
print(infiles)
for infile in infiles:
	if args.infile == "batch":
		prefix = infile.split(".")[0]
	else:
		prefix = args.outfile
	o = 0
	items = tblread(infile)
	#print(items)
	if len(items) > 1:
		print("Warning, more than one sequence detected")
		count = True
	else:
		count = False
	for i in items: 
		o += 1
		if count == True:
			outfile =  f"{prefix}_{o}.gff3"
		else:
			outfile =  f"{prefix}.gff3"
		with open(outfile, "w") as result:
			c = 0
			for a in items[i]:
				modules = items[i][a]
				#print(modules)
				contig = prefix
				if contig not in ("Sequence", "Name", "--------"): #this makes no sense
					start = modules["Range"][0]
					stop = modules["Range"][1]
					gene = a.replace(" ", "_")
					#oneletter = codechange[tRNA]
					score = 0
					#no introns assumed
					#include description!
					if modules["Strand"] == "forward":
						result.write("{0}\tMFannot\tgene\t{1}\t{2}\t{3}\t+\t.\tID={4}_{6};Name={4}_gene\n\
{0}\tMFannot\tCDS\t{1}\t{2}\t{3}\t+\t.\tID={4}{6};Name={4}\n\
{0}\tMFannot\texon\t{1}\t{2}\t{3}\t+\t.\tID={4}{6}_exon\n"
.format(contig, start, stop, score, gene, "", c))
						c += 1
					else:
						result.write("{0}\tMFannot\tgene\t{1}\t{2}\t{3}\t-\t.\tID={4}_{6};Name={4}_gene\n\
{0}\tMFannot\tCDS\t{1}\t{2}\t{3}\t-\t.\tID={4}_{6};Name={4}\n\
{0}\tMFannot\texon\t{1}\t{2}\t{3}\t-\t.\tID={4}_{6}_exon\n"
.format(contig, start, stop, score, gene, "", c))
						c += 1

	print("Converting done.")

