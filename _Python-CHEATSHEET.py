if isinstance(function, str): #determine if instance is string/list/set/whatever
print(os.path.isdir("/home/el"))
print(os.path.exists("/home/el/myfile.txt"))

#non-ascii character error:
import io
notincluded = io.open("abstracts.txt", "r", encoding="utf-8") #this is an analog of open(file).read()

import csv
csv.reader(open("cafe_table.tsv"), delimiter="\t", skipinitialspace=True)

break vs continue #break ends the for loop, continue ends the cycle for the current item
pass #use with if statements

import re
IPSpattern = r'"InterPro: (\w+)'
IPSID = re.search(IPSpattern, line).group(1)
for hit in re.findall(IPSpattern, i):
	blabla
#POZOR, re.match hleda na zacatku stringu

remove() vs discard()

import os
#to run os commands
if os.path.isfile("{}.nhr".format(args.database)) == False:
	os.system('makeblastdb -in {} -out {} -dbtype nucl -parse_seqids'.format(args.database, args.database))
else:
	print "skipping makeblastdb... database exists"

cmd = 'blastn'
db = 'VBRAnt'
query = 'CryptoDB-34_VbrassicaformisCCMP3155_AnnotatedProteins.fasta'
out = 'vbra_all_blast_mmetsp.xml'
evalue = 10
outfmt = 5
word_size = 4
threads = 4
maxevalue = 0.01
maxhits = ''
if os.path.isfile("vbra_all_blast_mmetsp.xml") == False:
	print('{} -query {} -db {} -out {} -evalue {} -outfmt {} -word_size {} -num_threads {}'.format(cmd, query, db, out, evalue, outfmt, word_size, threads))
	os.system('{} -query {} -db {} -out {} -evalue {} -outfmt {} -word_size {} -num_threads {}'.format(cmd, query, db, out, evalue, outfmt, word_size, threads))

#or
os.system('targetp -P {}-out.fasta > {}-targetp.txt'.format(prefix, prefix))

os.chdir("some/dir")
files = os.listdir('.')
files = [f for f in files if f.endswith(".fasta")]

#home = "/Users/zoliq/ownCloud/"
home = "/Volumes/zoliq data/OwnCloud/"
wd = home + "genomes/chromera/plastid proteome/"
os.chdir(wd)

if os.stat(fname).st_size > 0: #check for file size
	something()

directory_walker = os.walk('.')
for d in directory:
	#d[0] current directory string, d[1] subdirectories list, d[2] current directory file list
	print(d)

x = "NNNcagat"
prve_N = x.find('N')
posledne_N = x.rfind('N')
print(prve_N, posledne_N) #indexy


import argparse
#to parse arguments listed in the command after the script
parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-i', '--infile', help='Fasta/Phylip set to be trimmed', required=True)
parser.add_argument('-t', '--tree', help='Treefile for trimming', required=True)
parser.add_argument('-c', '--colour', help='Branch colours', default='all')
parser.add_argument('-n', '--num_seqs', help='maximum number of hits to recover', default='10')
parser.add_argument('-b', '--bootstrap', help='Boostrap calculation', action='store_true') # to only look for presence/absence of this argument


args = parser.parse_args()

infile = args.infile
num_seqs = args.num_seqs

##################################
#### Create working directory ####
##################################

if os.path.isdir("/Users/morpholino/OwnCloud/"):
	home = "/Users/morpholino/OwnCloud/"
elif os.path.isdir("/Volumes/zoliq data/OwnCloud/"):
	home = "/Volumes/zoliq data/OwnCloud/"
else:
	print("Please set a homedir")
if args.directory == ".":
	print("changing to default directory")
	defdir = "Robolab/phylo/"
	wd = home + defdir
	os.chdir(wd)
else:
	os.chdir(args.directory)


import sys
#communicating with the system
#another way of importing arguments from command line
with open(sys.argv[1]) as infile, open(sys.argv[2]) as outfile:
	something
#determine if operation system is OSX or linux
if sys.platform in ["darwin", "win32"]:
	choice = input().lower()
elif sys.platform.startswith("linux"):
	choice = raw_input().lower()
else:
	print("unrecognized OS, check how to use raw input")
	choice = raw_input().lower()

import time
print("Start time:", time.ctime())
print("Finish time:", time.ctime())



from collections import OrderedDict
#creates an ordered dictionary
seq_dict = OrderedDict()

#to avoid dictionary KeyErrors:
hightaxon = high_taxon_assignment_d.get(genus, "unassigned")


#read csv/tsv as dictionary:
import csv
with open('enz2annot.txt', 'r') as f:
	reader = csv.reader(f, delimiter='\t')
	enz2annot = {r[0]: r[1] for r in reader}


from Bio import SeqIO,AlignIO
#handling Fasta inputs
inFile = SeqIO.parse('input.txt', 'fasta')
from Bio.SeqUtils import GC
GC(rec.seq) #to calculate GC % of a sequence

with open('Split_fasta.txt', 'a') as result:
    for sequence in inFile:
        name = sequence.name
        seq = sequence.seq
        result.write('{}\t{}\n'.format(name,seq))
with open("safe-seq.fasta") as infile, open("safe-seq.phy", "w") as outfile:
	alignmentfile = AlignIO.read(infile, "fasta") 
	#apparently also works as follows:
	# alignmentfile = AlignIO.read("some_file.fasta", "fasta")
	AlignIO.write(alignmentfile, outfile, "phylip-relaxed")
#for batch data:
files = os.listdir('.')
files = [f.replace(".fasta","") for f in files if f.endswith(".fasta")]
for file in files:
	with open("{}.fasta".format(file)) as infile, open("{}.phy".format(file), "w") as outfile:
		alignmentfile = AlignIO.read(infile, "fasta")
		AlignIO.write(alignmentfile, outfile, "phylip-relaxed")
alignmentfile.get_alignment_length() #gets alignment length
alignmentfile[:, x] #get column of this index, watch out this is 0-based

with open('a', 'w') as a, open('b', 'w') as b: # open multiple files at the same time


from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
#handling sequences
mySL = coding_dna = Seq(args.SL, IUPAC.ambiguous_dna)
revSL = mySL.reverse_complement()
protein = mySL.translate()



#to remove bad chars from sequence names
#and other useful things
badchars = (",:.;()'")
fullname = ''.join(c for c in fullname if c not in badchars)
#OR
prediction = re.sub('[():,]', '', prediction)


if sequence.startswith('>'):
	counter += 1


#to print shortened decimals
print('%.3f' % (somenumber))
print("{:.2f}".format(somenumber))


#avoid trailing newline in fastq
for index, item in enumerate(reads):
	if index != len(reads) - 1:
		name = item[0].split(' ')[0]
		out.write('{}/{}\n{}\n'.format(name, direction, '\n'.join(item[1:])))
	else:
		name = item[0].split(' ')[0]
		out.write('{}/{}\n{}\n+\n{}'.format(name, direction, item[1], item[3]))	


from tqdm import tqdm
for f in tqdm(files, desc="progress"):
	do whatever scripting

#or
with open(queryfile) as f:
	seqcount = f.read().count(">")
print("Sequence reading progress: {:.1f}%".format(100*(c / seqcount)))

############################
####  useful functions  ####
############################

#calculate centroid of a protein
from __future__ import division
import numpy as np

data = np.genfromtxt('file.pdb') #effective way how to import tables directly into an array
array = data[:, -3:] #if the data is here
np.mean(data[:,-3:], axis=0) #mean along the vertical axis, more or less means calculate centroid


def query_yes_no(question, default="yes"):
	"""Ask a yes/no question via raw_input() and return their answer.

	"question" is a string that is presented to the user.
	"default" is the presumed answer if the user just hits <Enter>.
		It must be "yes" (the default), "no" or None (meaning
		an answer is required of the user).

	The "answer" return value is True for "yes" or False for "no".
	"""
	valid = {"yes": True, "y": True, "ye": True,
			 "no": False, "n": False}
	if default is None:
		prompt = " [y/n] "
	elif default == "yes":
		prompt = " [Y/n] "
	elif default == "no":
		prompt = " [y/N] "
	else:
		raise ValueError("invalid default answer: '{}'" .format(default))

	while True:
		sys.stdout.write(question + prompt)
		if sys.platform in ["darwin", "win32"]:
			choice = input().lower()
		elif sys.platform.startswith("linux"):
			choice = raw_input().lower()
		else:
			print("unrecognized OS, check how to use raw input")
			choice = raw_input().lower()

		if default is not None and choice == '':
			return valid[default]
		elif choice in valid:
			return valid[choice]
		else:
			sys.stdout.write("Please respond with 'yes' or 'no' "
							 "(or 'y' or 'n').\n")

def delbadchars(string):
	""" remove unneeded characters """
	badchars = ("|+,:;()' ") #also []/@
	n = []
	for c in string:
		if c in badchars:
			c = "_"
		n.append(c)
	result = "".join(n)
	return result


def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""

    # from whichcraft import which
    from shutil import which

    return which(name) is not None


def overlaps(moduleA, moduleB):
	startA = moduleA[0]
	endA = moduleA[1]
	startB = moduleB[0]
	endB = moduleB[1]
	if startA <= startB <= endA:
		if startA <= endB <= endA:
			overlap = "embed"
		else:
			overlap = "overlap"
	elif startA <= endB <= endA:
		overlap = "overlap"
	elif startB <= startA <= endB and startB <= startB <= endB:
		overlap = "embed"
	else:
		overlap = "none"
	return overlap


def decor(string):
	def wrap():
		print("===============")
		print(string)
		print("===============")
	return wrap