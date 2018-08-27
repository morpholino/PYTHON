import os
from Bio import SeqIO
from Bio.Seq import Seq

homedir = "/Users/zoliq/Downloads/"
#homedir = "/Volumes/zoliq data/ownCloud"
wd = homedir + "."
os.chdir(wd)

infasta = SeqIO.parse("elonga_allreads.fasta", "fasta")
counter = 0
baits = ("CTGCTGGACGCTGTCCTGCTTTGAGATGAAGATGTGGTGGGTGATCCTGTGGACGTCCTACTGCAGGCGGTGGTCCGGGCGCAGCTGGCAGTGGCCGCCTG",
		"TCGGCTGCCACCCGACGGCCCAGAAGCCCCCCAGCTGGATGAAGCGCTGCGCCGGCTCCTACAGGGACACTATAGCCACGTCGCCTTCACCTCGCGGGCGG",
		"ACCCGACGGCCCAGAAGCCCCCCAGCTGGATGAAGCGCTGCGCCGGCTCCTACAGGGACACTATAGCCACGTCGCCTTCACCTCGCGGGCGGGCATCGAGG", 
		"CTATAGCCACGTCGCCTTCCCCTCGCGGGCGGGCATCGAGGCGGTGCTGGAGCGGCTGGAGGCTCTGGGGCGGGGTGAACCGCCGTCCCGATTCCCCGTTG",
		"CGGCTGGAGGCTCTGGGGCGGGGTGAACCGCCGTCCCGATTCCCCGTTGGTCAGGATCGGCCCGCAGCACTCGCAGCAGGTGGGGCGTTGTTGGCCCAAAC",
		"GGGGCGGGGTGAACCGCCGTCCCGATTCCCCGTTGGTCAGGATCGGCCCGCAGCACTCGCAGCAGGTGGGGCGTTGTTGGCCCAAACCAACGTCAGTGTGT")
# first bait is for control if the script works, it is the fifth read in the dataset
for seq in infasta:
	counter += 1
	if counter % 1000000 == 0:
		print(counter)
	# if seq.seq in baits:
	# 	print(">{}\n{}\n".format(seq.description, seq.seq))
	# elif seq.seq.reverse_complement() in baits:
	# 	print(">{}\n{}\n".format(seq.description, seq.seq))

	if seq.name in ("HWI-ST143:591:C16F2ACXX:4:1203:2666:76392",
		"HWI-ST143:591:C16F2ACXX:4:2111:13138:27359", 
		"HWI-ST143:591:C16F2ACXX:4:1103:19968:31125",
		"HWI-ST143:591:C16F2ACXX:4:2210:9576:68454",
		"HWI-ST143:591:C16F2ACXX:4:2308:15402:37115"):
		print(">{}\n{}\n".format(seq.description, seq.seq))