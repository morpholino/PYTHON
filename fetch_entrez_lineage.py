from ete3 import NCBITaxa
#http://etetoolkit.org/docs/2.3/tutorial/tutorial_ncbitaxonomy.html
ncbi = NCBITaxa()

"""with open("fetch_entrez_missing.txt") as f:
	missing = f.read().split("\n")
#print(missing)
name2taxid = ncbi.get_name_translator(missing)
print(name2taxid)
#TESTING PURPOSES ONLY:
#missing = ['Homo', 'Aspergillus', 'Haloquadratum']
#taxid2name = ncbi.get_taxid_translator([9606, 9443])
#rankofnode = ncbi.get_rank([9606, 9443])

for genus in name2taxid:
	#retriev at least one species:
	descendants = ncbi.get_descendant_taxa(genus)
	lineage = ncbi.get_lineage(descendants[0])[2:7]
	names = ncbi.get_taxid_translator(lineage)
	rank = [names[taxid] for taxid in lineage]
	if "Eukaryota" in rank:
		rank.remove("Eukaryota")
	print("{}\t{}\t{}".format(genus, ";".join(rank)))"""

with open("/Users/morpholino/OwnCloud/AndyLab/MMETSP_mcl/species_list.txt") as f:
	missing = f.read().split("\n")
	missing = [x.split("\t")[1] for x in missing] #this contains species taxids now
print(missing)
name2taxid = ncbi.get_name_translator(missing)
print(name2taxid)

for genus in missing: #name2taxid:
	#retriev at least one species:
	descendants = ncbi.get_descendant_taxa(genus)
	lineage = ncbi.get_lineage(descendants[0])[2:] #[2:7]
	names = ncbi.get_taxid_translator(lineage)
	rank = [names[taxid] for taxid in lineage]
	#if "Eukaryota" in rank:
	#	rank.remove("Eukaryota")
	print("{}\t{}".format(genus, ";".join(rank))) #name2taxid[genus][0], 

