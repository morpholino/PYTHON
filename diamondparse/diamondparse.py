from Bio import SeqIO
import os
import argparse

if not os.path.isdir("diamond_out"):
	os.system("mkdir diamond_out")
	print("DIAMONDPARSE: Creating target directory")

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-p', '--pool', help='Pool results?', default=True)

args = parser.parse_args()


print("DIAMONDPARSE: Starting...")
#FILENAMES:
BLASTOUT = 'diamond.out'
SORTDIR = 'diamond_out/'
QUERYFILE = 'diamondquery.fasta'
LOG = 'Diamond-mod-log.txt' #NOTE this is for a different than default-setting diamond outfiles
#evaluethresh = 0.0001 #e-value not in output

seq_count = 0
sequence_dict = {}
for sequence in SeqIO.parse(QUERYFILE, 'fasta'):
	#print(sequence.name)
	seq_count += 1
	sequence_dict[sequence.name] = {} #{sequence.description: sequence.seq} if query needs to be included
print("DIAMONDPARSE: Sequences extracted from query file: {}".format(seq_count))

#goodhits = set()
#badqueries = set()
#with open('example-BLAST.txt') as infile:
with open(BLASTOUT) as infile:
	for line in infile:
		line = line.strip().split("\t")
		query = line[0]
		if query not in sequence_dict:
			print("ERROR! unspecified query:", query)
		else:
			sequence_dict[query].update({line[1]: (line[2], line[-1])})

if args.pool in {True, "True", "true"}:
	print("DIAMONDPARSE: pooling hit sequences to {}{}...".format(SORTDIR, "pool"))
	with open("{}{}_out.fasta".format(SORTDIR, "pool"), "w") as result: 
		for query in sequence_dict:
			for item in sequence_dict[query]:
				#print(item)
				hit = sequence_dict[query][item]
				#print(hitseq)
				#print(item, hit)
				result.write(">{}\n{}\n".format(hit[0], hit[1]))

else:
	print("DIAMONDPARSE: pooling disabled, hit sequences in {}".format(SORTDIR))
	for query in sequence_dict:
		with open("{}{}_out.fasta".format(SORTDIR, query), "w") as result: 
			for item in sequence_dict[query]:
				#print(item)
				hit = sequence_dict[query][item]
				#print(hitseq)
				#print(item, hit)
				result.write(">{}\n{}\n".format(hit[0], hit[1]))

print("DIAMONDPARSE: Extracted proteins written to files. Finished.")
