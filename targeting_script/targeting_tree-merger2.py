#function to read fasta and create a list of sequences
print("This script accepts prediction outputs, along with a phylobayes tree.")
print("The script then outputs the renamed treefile with full taxa names and the prediction added at the end.")
print("Please, enter the prefix (-p) of prediction files to be analyzed...(xxx.txt)")
print("Treefiles are read automatically from the current folder, or can be entered with (-t).")
print("Optionally, provide a taxa replacement key file in a .tsv format (-a)\n###############################################")

import argparse
import os
import re
from Bio import SeqIO

#### Change to workdir ####
###########################

homedir = "/Users/zoliq/ownCloud/"
#homedir = "/Volumes/zoliq data/ownCloud"
wd = homedir + "genomes/phatr/phatr mitoglyco/huge alignments/PASTA alignments nonconverging/good sets/localization/"
os.chdir(wd)

#### Collect Input ####
#######################

parser = argparse.ArgumentParser(description='How to use argparse')
parser.add_argument('-p', '--prefix', help='Prediction files prefix', default='seq_eukaryotes')
parser.add_argument('-t', '--treefile', help='Treefile', default='none')
parser.add_argument('-a', '--accessions', help='Accession rename key file', default='none')

args = parser.parse_args()

prefix = args.prefix
treefile = args.treefile
accessions = args.accessions
#ONLY FOR TEST PURPOSES:
accessions = "leaf_renaming.txt"

print("FILE prefix defined: %s" %(prefix))
if accessions != 'predefined':
	print("Taxa replacement key defined: %s" %(accessions))
else:
	print("Only predefined codes are used.")

#### Functions ####
###################

def second_largest(numbers):
    count = 0
    m1 = m2 = float('-inf')
    for x in numbers:
        count += 1
        if x > m2:
            if x >= m1:
                m1, m2 = x, m1            
            else:
                m2 = x
    return m2 if count >= 2 else None


#### Open and parse inputs ####
###############################

preds_d = {} #predictions dictionary
#preds_d['accession'] = {'hectar':'ND', 'signalp': 'ND', 'asafind': 'ND', 'targetp': 'ND', 'ML2ANIM': 'ND', 'ML2PLANT': 'ND'}


asafind = open(prefix + "-asafind.txt").read().split('\n')
possiblepredsasafind = {'No idea': 'N/A', 'Yes, Low': 'PLASTID, LOW', 'NA': 'NON-PLASTID', 'No': 'NON-PLASTID', 'Yes, High': 'PLASTID, HIGH'}
for item in asafind:
	#Identifier	SignalP	ASAfind cleavage position	ASAfind/SignalP cleavage site offset	ASAfind 20aa transit score	ASAfind Prediction
	#[0]		[1]		[2]							[3]										[4]							[5]	
	
	if not item.startswith('#') and len(item) != 0:
		item = item.split('\t')
		name = item[0]
		pred = item[5]
		pred = possiblepredsasafind[pred]
		preds_d[name] = {'asafind': pred}


hectar = open(prefix + "-hectar.txt").read().split('\n')
possiblepredshectar = {'mitochondrion': 'MITOCHONDRION', 'no signal peptide or anchor': 'OTHER', 'other localisation': 'OTHER', 'chloroplast': 'PLASTID', 'signal anchor': 'SIGNAL', 'signal peptide': 'SIGNAL'}
for item in hectar:
	#protein id	predicted targeting category	signal peptide score	type II signal anchor score	chloroplast score	mitochondrion score	other score
	#[0]		[1]								[2]						[3]							[4]					[5]					[6]
	#BUT!
	#protein id	step	error message	signal peptide score	type II signal anchor score
	#[0]		[1]		[2]				[3]						[4]

	if not item.startswith('#') and len(item) != 0:
		item = item.split('\t')
		item = list(filter(None, item)) #remove list() function in python2.x
		name = item[0]
		pred = item[1]

		if pred not in ("Failed(U/X)", "error message"):
			#STRAMENOPILES PREDICTIONS:
			#protein id	predicted targeting category	signal peptide score	type II signal anchor score		chloroplast score	mitochondrion score		other score
			#[0]		[1]								[2]						[3]								[4]					[5]						[6]
			
			#OPISTHOKONTS PREDICTIONS:
			#protein id	predicted targeting category	signal peptide score	type II signal anchor score		mitochondrion score		other score
			#[0]		[1]								[2]						[3]								[4]						[5]					

			#UNSPECIFIED EUKARYOTES PREDICTIONS:
			#protein id	predicted targeting category	signal peptide score	type II signal anchor score		other score
			#[0]		[1]								[2]						[3]								[4]
			
			#ERROR:
			#protein id	step	error message	signal peptide score	type II signal anchor score
			#[0]				[1]				[2]						[3]
			pred = possiblepredshectar[pred]
			if len(item) == 7:
				#HECTARPREDICTION "stramenopiles"
				item = [k.replace('-', '0') for k in item]
				hectarpreds_l = [float(item[2]), float(item[4]), float(item[5]), float(item[6])]
				hectarpreds_d = {float(item[2]): 'SP', float(item[4]): 'cTP', float(item[5]): 'mTP', float(item[6]): 'other'}
				pred = ("{}_({}:{} > {}:{})".format(pred, hectarpreds_d[max(hectarpreds_l)], max(hectarpreds_l), hectarpreds_d[second_largest(hectarpreds_l)], second_largest(hectarpreds_l)))
			elif len(item) == 6:
				#HECTARPREDICTION "animals+fungi"
				item = [k.replace('-', '0') for k in item]
				hectarpreds_l = [float(item[2]), float(item[4]), float(item[5])]
				hectarpreds_d = {float(item[2]): 'SP', float(item[4]): 'mTP', float(item[5]): 'other'}
				pred = ("{}_({}:{} > {}:{})".format(pred, hectarpreds_d[max(hectarpreds_l)], max(hectarpreds_l), hectarpreds_d[second_largest(hectarpreds_l)], second_largest(hectarpreds_l)))
			elif len(item) == 6:
				#HECTARPREDICTION "eukaryotes"
				item = [k.replace('-', '0') for k in item]
				hectarpreds_l = [float(item[2]), float(item[4])]
				hectarpreds_d = {float(item[2]): 'SP', float(item[4]): 'other'}
				pred = ("{}_({}:{} > {}:{})".format(pred, hectarpreds_d[max(hectarpreds_l)], max(hectarpreds_l), hectarpreds_d[second_largest(hectarpreds_l)], second_largest(hectarpreds_l)))
			curdict = {'hectar': pred}
			try:
				preds_d[name].update(curdict)
			except KeyError:
				print("Hectar: bad key " + name)
				preds_d[name] = {'hectar': pred}
		else:
			pred = "undefined amino acids U/X"
			try:
				preds_d[name].update(curdict)
			except KeyError:
				print("Hectar: bad key " + name)
				preds_d[name] = {'hectar': pred}


ML2ANIM = open(prefix + "-ML2animalHI.txt").read().split('\n')
possiblepredsml2 = {"chloroplast": "PLASTID", 'cytoplasmic': "CYTOSOL", "ER": "ER", "extracellular": "EXTRACELLULAR", "Golgi apparatus": "GOLGI", "mitochondrial": "MITOCHONDRION", "nuclear": "NUCLEUS", "lysosomal": "LYSOSOME", "peroxisomal": "PEROXISOME", "plasma membrane": "PLASMA MEMBRANE", "vacuolar": "VACUOLE"}
for item in ML2ANIM:
	#protein	decreasing localization predictions from: cytoplasmic, ER, extracellular, Golgi apparatus, mitochondrial, nuclear, lysosomal, peroxisomal, plasma membrane
	if not item.startswith('#') and len(item) != 0:
		item = item.split('\t')
		name = item[0]
		Loc = possiblepredsml2[item[1].split(':')[0]]
		pred = ("{}_({} > {})".format(Loc,item[1], item[2]))
		try:
			preds_d[name].update({"ML2ANIMAL": pred})
		except KeyError:
			print("MultiLoc2-animal: bad key " + name)
			preds_d[name] = {'ML2ANIMAL': pred}



ML2PLANT = open(prefix + "-ML2plantHI.txt").read().split('\n')
for item in ML2PLANT:
	#protein	decreasing localization predictions from: chloroplast, cytoplasmic, ER, extracellular, Golgi apparatus, mitochondrial, nuclear, lysosomal, peroxisomal, plasma membrane, vacuolar
	if not item.startswith('#') and len(item) != 0:
		item = item.split('\t')
		name = item[0]
		Loc = possiblepredsml2[item[1].split(':')[0]]
		pred = ("{}_({} > {})".format(Loc,item[1], item[2]))
		try:
			preds_d[name].update({"ML2PLANT": pred})
		except KeyError:
			print("MultiLoc2-plant: bad key " + name)
			preds_d[name] = {'ML2PLANT': pred}


signalp = open(prefix + "-signalp.txt").read().split('\n') #sensitive version 4.1
for item in signalp:
	# SignalP-3.0 euk predictions  
	# name 	Cmax 	pos ? 	Ymax	pos ? 	Smax 	pos ? 	Smean 	?	D 	? 	name 	!	Cmax	pos ? 	Sprob	?
	#  0	1		2	3	4		5	6	7		8	9	10		11	12	13	14		15	16		17	18	19		20
	# SignalP-4.1 euk predictions  
	# names Cmax	pos 	Ymax	pos 	Smax	pos 	Smean 	D 	?	Dmaxcut 	Networks-used 
	#  0	1		2		3		4		5		6		7		8	9	10			11

	if not item.startswith('#') and len(item) != 0:
		item = item.split()
		name = item[0]
		if len(item) == 21:
			#SIGNALPVERSION 3.0
			pred = item[13]
			pred = pred.replace("N", "OTHER")
			pred = pred.replace("Y", "SIGNAL")
		if len(item) == 12:
			#SIGNALPVERSION 4.x
			pred = item[9]
			pred = pred.replace("N", "OTHER")
			pred = pred.replace("Y", "SIGNAL")
		try:
			preds_d[name].update({'signalp': pred})
		except KeyError:
			print(item)
			print("SignalP: bad key " + name)
			preds_d[name] = {'signalp': pred}


targetp = open(prefix + "-targetp.txt").read().split('\n') #PLANT + NONPLANT
possiblepredstargetp = {"M": "MITOCHONDRION", "C": "PLASTID", "S": "SIGNAL", "_": "OTHER"}
for item in targetp:
	#Name       Len     mTP     SP      other   Loc     RC
	#[0]		[1]		[2]		[3]		[4]		[5]		[6]
	#BUT PLANT!
	#Name       Len     cTP     mTP     SP      other   Loc     RC
	#[0]		[1]		[2]		[3]		[4]		[5]		[6]		[7]

	if len(item.split()) > 1 and item.split()[1].isnumeric(): #takes only lines with prediction
		item = item.split()
		name = item[0]
		if len(item) == 7:
			targetpreds_l = [float(item[2]), float(item[3]), float(item[4])]
			targetpreds_d = {float(item[2]): 'mTP', float(item[3]): 'SP', float(item[4]): 'other'}
			Loc = possiblepredstargetp[item[5]]
			pred = ("{}_({}:{} > {}:{})".format(Loc, targetpreds_d[max(targetpreds_l)], max(targetpreds_l), targetpreds_d[second_largest(targetpreds_l)], second_largest(targetpreds_l)))
		elif len(item) == 8:
			targetpreds_l = [float(item[2]), float(item[3]), float(item[4]), float(item[5])]
			targetpreds_d = {float(item[2]): 'cTP', float(item[3]): 'mTP', float(item[4]): 'SP', float(item[5]): 'other'}
			Loc = possiblepredstargetp[item[6]]
			pred = ("{}_({}:{} > {}:{})".format(Loc, targetpreds_d[max(targetpreds_l)], max(targetpreds_l), targetpreds_d[second_largest(targetpreds_l)], second_largest(targetpreds_l)))
		try:
			preds_d[name].update({'targetp': pred})
		except KeyError:
			print("TargetP: bad key " + name)
			preds_d[name] = {'targetp': pred}
	else:
		pass

print("preds_dictionary collected")


#### Read own codes ####
########################
#read predefined taxa codes
taxa = {
	'crypGUILt': 'Guillardia theta', 'rhiAMMOs': 'Ammonia sp.', 'dinGYMNc': 'Gymnodinium catenatum', 'funASPEf': 'Aspergillus fumigatus', 
	'kytORYZs': 'Oryza sativa Japonica', 'amoPAULc': 'Paulinella chromatophora', 'dinLINGp': 'Lingulodinium polyedra', 'oomPHYra': 'Phytophtora ramorum', 
	'hwal': 'Haloquadratum walsbyi', 'alfaNITRh': 'Nitrobacter hamburgensis', 'chryHETak': 'Heterosigma akashiwo', 'mlep': 'Mycobacterium leprae', 
	'chlaLOTam': 'Lotharella amoebiformis', 'strTHRAs': 'Thraustochytrium sp.', 'rhiPLASb': 'Plasmodiophora brassicae', 
	'strNANNO': 'Nannochloropsis gaditana B-31', 'redTIMSo': 'Timspurckia oligopyrenoides', 'chryCHATs': 'Chattonella subsalsa', 
	'rhiROSAs': 'Rosalina sp.', 'tgon': 'Toxoplasma gondii', 'ddis': 'Dictyostelium discoideum', 'rhiSORIs': 'Sorites sp.', 
	'strCAFca': 'Cafeteria Caron', 'mmar': 'Methanococcus maripaludi', 'perkPmari': 'Perkinsus marinus', 'dhan': 'Debaryomyces hansenii', 
	'crei': 'Chlamydomonas reinhardtii', 'hsap': 'Homo sapiens', 'gamaYERSp': 'Yersinia pestis', 'sent': 'Salmonella enterica', 
	'firmLISTm': 'Listeria monocytogenes', 'bant': 'Bacillus anthracis', 'osat': 'Oryza sativa Japonica', 'grnBATHp': 'Bathycoccus prasinos', 
	'apiCRYPm': 'Cryptosporidium muris', 'redGALDs': 'Galdieria sulphuraria', 'actiMYCOt': 'Mycobacterium tuberculosis', 'dinHETtr': 'Heterocapsa triquestra', 
	'cele': 'Caenorhabditis elegans', 'cbot': 'Clostridium botulinum', 'rhiPARTg': 'Partenskyella glossopodia', 'dinCRYPc': 'Crypthecodinium cohnii', 
	'cilFAVta': 'Favella taraikaensis', 'deth': 'Dehalococcoides ethenogenes', 'atha': 'Arabidopsis thaliana', 'cilTETRt': 'Tetrahymena thermophila', 
	'cjej': 'Campylobacter jejuni', 'ypes': 'Yersinia pestis', 'archPYROa': 'Pyrobaculum aerophilum', 'redRHSOm': 'Rhodosorus marinus', 
	'redCHONc': 'Chondrus crispus', 'archSULFt': 'Sulfolobus tokodaii', 'cneg': 'Cryptococcus neoformans', 'strECTsi': 'Ectocarpus siliculosus', 
	'cyanSYNEs': 'Synechococcus sp.', 'alfaRICKc': 'Rickettsia conori', 'grnPOLYp': 'Polytomella parva', 'wend': 'Wolbachia endosymbiont', 
	'gamaSHEWb': 'Shewanella baltica', 'oomPYTHu': 'Pythium ultimum', 'chryPTERd': 'Pteridomonas danica', 'aaeo': 'Aquifex aeolicus', 
	'cyanPROCm': 'Prochlorococcus marinus', 'tpal': 'Treponema pallidum', 'redCYANm': 'Cyanidioschyzon merolae', 'haptPRYMp': 'Prymnesium parvum', 
	'betaCUPRn': 'Cupriavidus necator', 'dinCERAf': 'Ceratium fusus', 'cpne': 'Chlamydophila pneumoniae', 'haptPLEUc': 'Pleurochrysis carterae', 
	'strAURAl': 'Aurantiochytrium limacinum', 'glauGLOwi': 'Gloeochaete wittrockiana', 'apiGREGn': 'Gregarina niphandrodes', 
	'dinAMPHm': 'Amphidinium massartii', 'cilEUPfo': 'Euplotes focardii', 'chlaLOTgl': 'Lotharella globosa', 'dinAMPHc': 'Amphidinium carterae', 
	'rhiMINch': 'Minchinia chitonis', 'grnPYRpa': 'Pyramimonas parkeae', 'tbru': 'Trypanosoma brucei', 'vcho': 'Vibrio cholerae', 
	'cyanTHERe': 'Thermosynechococcus elongatus', 'gsul': 'Geobacter sulfurreducens', 'dmel': 'Drosophila melanogaster', 'redPORpu': 'Porphyra purpurea', 
	'grnCHLAr': 'Chlamydomonas reinhardtii', 'cyanCYANp': 'Cyanothece sp.', 'dinOXYma': 'Oxyrrhis marina', 
	'funLACCb': 'Laccaria bicolor', 'excNAEgr': 'Naegleria gruberi', 'tcru': 'Trypanosoma cruzi', 'strTHAps': 'Thalassiosira pseudonana', 
	'mtub': 'Mycobacterium tuberculosis', 'dinALEXc': 'Alexandrium catenella', 'excPERco': 'Percolomonas cosmopolitus', 'dinALEXa': 'Alexandrium andersonii', 
	'cilLITOp': 'Litonotus pictus', 'alfaMAGNm': 'Magnetospirillum magneticum', 'grnCHLOv': 'Chlorella variabilis', 'cpos': 'Coccidioides posadasii', 
	'strAPLAs': 'Aplanochytrium stocchinoi', 'redRHOma': 'Rhodella maculata', 'grnOSTRt': 'Ostreococcus tauri', 'eugEUTn': 'Eutreptiella gymnastica NIES-381', 
	'grnVOLVc': 'Volvox carteri f. nagariensis', 'aniMONOb': 'Monosiga brevicollis', 'chryOCHRO': 'Ochromonas sp.', 'sfle': 'Shigella flexneri', 
	'eugEUTc': 'Eutreptiella gymnastica CCMP1594', 'firmSTAPa': 'Staphyllococcus aureus', 'dinKARmi': 'Karlodinium micrum', 'sman': 'Schistosoma mansoni', 
	'chlaLOToc': 'Lotharella oceanica', 'cilPARAt': 'Paramecium tetraurelia', 'dinDURIb': 'Durinskia baltica', 'apiCHROv': 'Chromera velia', 
	'kytPHYSp': 'Physcomitrella patens', 'grnCOCCs': 'Coccomyxa subellipsoidea', 'cilFABRs': 'Fabrea salina', 'spom': 'Schizosaccharomyces pombe', 
	'ylip': 'Yarrowia lipolytica', 'dinPYROb': 'Pyrodinium bahamense', 'strCAFro': 'Cafeteria roenbergensis', 'amoACANc': 'Acanthamoeba castellanii', 
	'kytARABt': 'Arabidopsis thaliana', 'dinKARb': 'Karenia brevis', 'ggal': 'Gallus gallus', 'betaRALSs': 'Ralstonia solanacearum', 
	'crypRHOsa': 'Rhodomonas salina', 'lmon': 'Listeria monocytogenes', 'apiVOROp': 'Voromonas pontica', 
	'actiCORYd': 'Corynebacter diphteriae', 'betaVERMe': 'Verminephrobacter eiseniae', 'rhiELPHm': 'Elphidium margitaceum', 'ctep': 'Chlorobium tepidum', 
	'alfaAZOSs': 'Azospirillum sp.', 'cilTIARf': 'Tiarina fusus', 'chlaBIGna': 'Bigelowiella natans', 'eugEUGlo': 'Euglena longa', 'ppat': 'Physcomitrella patens', 
	'dinNOCTs': 'Noctiluca scintillans', 'haptEMILh': 'Emiliania huxleyi', 'cyanLYNGp': 'Lyngbya sp.', 'cyanANABv': 'Anabaena variabilis', 
	'calb': 'Candida albicans', 'apiTOXOg': 'Toxoplasma gondii', 'cper': 'Clostridium perfringens', 'tmar': 'Thermotoga maritima', 'wsuc': 'Wolinella succinogenes', 
	'pfal': 'Plasmodium falciparum', 'ecol': 'Escherichia coli', 'chlaLOTHs': 'Lotharella sp.', 'betaBURKc': 'Burkholderia cenocepacia', 
	'cilPLATm': 'Platyophrya macrostoma', 'eugEUGgr': 'Euglena gracilis', 'redPORae': 'Porphyridium aerugineum', 
	'bcidFLAVc': 'Flavobacterium columnare', 'mmus': 'Mus musculus', 'glauCYPTg': 'Cyanoptyche gloeocystis', 'mjan': 'Methanocaldococcus jannaschii', 
	'grnPYRam': 'Pyramimonas amylifera', 'ssol': 'Sulfolobus solfataricus', 'cyanNODUs': 'Nodularia spumigena', 'grnDUNte': 'Dunaliella tertiolecta', 
	'actiSTREc': 'Streptomyces coelicolor', 'glauCYApa': 'Cyanophora paradoxa', 'chryVAUCH': 'Vaucheria litorea', 'grnMICRp': 'Micromonas pusilla', 
	'bcidBACTf': 'Bacteroides fragilis', 'apiTHEIe': 'Theileria equi', 'eugTRYPb': 'Trypanosoma brucei', 'aory': 'Aspergillus oryzae', 
	'strPHAtr': 'Phaeodactylum tricornutum', 'dinKRYPf': 'Kryptoperidinium foliaceum', 'amoDICTd': 'Dictyostelium discoideum', 'firmBACIa': 'Bacillus anthracis', 
	'saur': 'Staphylococcus aureus', 'dinHETro': 'Heterocapsa rotundata', 'msed': 'Metallosphaera sedula', 'dinGLENf': 'Glenodinium foliaceum', 
	'tvol': 'Thermoplasma volcanium', 'archTHERv': 'Thermoplasma volcanium', 'apiVITbr': 'Vitrella brassicaformis', 'spne': 'Streptococcus pneumoniae', 
	'bcidPORPg': 'Porphyromonas gingivalis', 'redERYTm': 'Erythrolobus madagascarensis', 'bmaa': 'Brugia malayi', 'bmal': 'Burkholderia mallei', 
	'cilCONDm': 'Condylostoma magnum', 'alfaRHODs': 'Rhodobacter sphaeroides', 'eugLEISm': 'Leishmania major', 'ncra': 'Neurospora crassa', 
	'bsui': 'Brucella suis', 'crypGUIth': 'Guillardia theta', 'mmul': 'Macaca mulatta', 'kytSELAm': 'Selaginella moellendorffii', 'dinSYMBs': 'Symbiodinium sp.', 
	'gamaVIBRc': 'Vibrio cholerae', 'rsol': 'Ralstonia solanacearum', 'apiEIMEt': 'Eimeria tenella', 'strSKEma': 'Skeletonema marinoi', 
	'aful': 'Archaeoglobus fulgidus', 'afum': 'Aspergillus fumigatus', 'cbur': 'Coxiella burnetii', 'archPICRt': 'Picrophilus torridus', 
	'strAUREa': 'Aureococcus anophagefferens', 'drad': 'Deinococcus radiodurans', 'egos': 'Eremothecium gossypii', 'klac': 'Kluyveromyces lactis', 
	'cyanCROCw': 'Crocosphaera watsonii', 'ftul': 'Francisella tularensis', 'rbal': 'Rhodopirellula baltica', 'redMADAe': 'Madagascaria erythrocladoides', 
	'apiNEOSc': 'Neospora caninum', 'amoENTAh': 'Entamoeba histolytica', 'crypCRYpa': 'Cryptomonas paramecium', 'ecab': 'Equus caballus', 
	'funCRYPn': 'Cryptococcus neoformans', 'dinPRORm': 'Prorocentrum minimum', 'cimm': 'Coccidioides immitis', 'scer': 'Saccharomyces cerevisiae', 
	'chryDINOB': 'Dinobryon sp.', 'redERYTa': 'Erythrolobus australicus', 'aniDROSm': 'Drosophila melanogaster', 'bcidPREVr': 'Prevotella ruminicola', 
	'strNANNg': 'Nannochloropsis gaditana', 'dinDINac': 'Dinophysis acuminata'
}
#addition of specific taxacodes entered by -a to the taxa vocabulary
if accessions != 'none':
	codes = open(accessions).read().split("\n")
	#codes = open("leaf_renaming.txt").read().split("\n")
	for code in codes:
		if len(code) != 0:
			code = code.split("\t")
			taxa[code[0]] = code[1]

#### Leafsearch ####
####################
"""
#UNCOMMENT THIS PART TO PREPARE A TAXA-REPLACEMENT FILE
leaves = open("leaf_names.txt").read().split("\n")

finalnameset = set()
nohigherordername = {}
inFasta = SeqIO.parse("allseq.fasta", 'fasta')
for seq in inFasta:
	finalnameset.add(seq.description)
	nohigherordername[seq.description.split("@")[0]] = seq.description

leafrenaming = {}
count = 0
for key in preds_d:
	modkey = key.split("@")[0]
	nohigherordername[modkey] = key
for leaf in leaves:
	if len(leaf) != 0:
		origleafname = leaf.split()[1]
		finalleafname = leaf.split()[1].replace("&"," ")
		if finalleafname in finalnameset:
			leafrenaming[origleafname] = finalleafname
		elif finalleafname.split("@")[0] in nohigherordername:
			leafrenaming[origleafname] = finalleafname
		predleafname = leaf.split()[1].replace("&","_")
		if predleafname in preds_d:
			leafrenaming[origleafname] = predleafname
		elif predleafname.split("@")[0] in nohigherordername:
			leafrenaming[origleafname] = nohigherordername[predleafname.split("@")[0]]
for key in preds_d:
	if key not in leafrenaming.values():
		print(key)
		count += 1


with open("another_leaf_renaming.txt", "w") as result:
	for origname in leafrenaming:
		result.write("{}\t{}\n".format(origname, leafrenaming[origname]))
print("To find manually: ", count)

count = 0
leaf_renaming = open("leaf_renaming.txt").read().split("\n")
for line in leaf_renaming:
	if len(line) != 0:
		line = line.split("\t")
		taxa[line[0]] = line[1]
for leaf in leaves:
	if len(leaf) != 0:
		leaf = leaf.split()[1]
		if leaf not in taxa:
			print(leaf)
			count += 1

print("To find manually: ", count)
"""


#### Main ####
##############

heterotrophs = {'Ciliophora', 'Stereomyxa', 'Tubulinea', 'Metazoa', 'Hilomonadea', 'Lobosa', 'Metamonada', 'Fungi', 
'Choanoflagellida', 'Foraminifera', 'Conosa', 'Amoebozoa', 'Excavata'}
opisthokonts = {'Metazoa', 'Fungi', 'Choanoflagellida'} 	#Hectar / MultiLoc-animal
otherhetero = heterotrophs - opisthokonts 					#TargetP / MultiLoc-animal
primary = {'Chlorophyta', 'Rhodophyta', 'Streptophyta', 'Glaucocystophyta', 'Glaucophyta'}	#TargetP / MultiLoc-plant
higherorder = {'Dinophyta', 'Discoba', 'Haptophyta', 'Picobiliphyta', 'Synurophyceae', 'Alveolata_X', 'Cryptophyta', 
'Apicomplexa', 'Cercozoa', 'Alveolata', 'Hacrobia_X', 'Perkinsea', 'Stramenopiles_X', 'Chromerida'}
eukaryote = heterotrophs | primary | higherorder
stramenopiles = {'Synurophyceae', 'Stramenopiles_X'}		#Hectar / ASAFind
otherhigher = higherorder - stramenopiles 					#TargetP / ASAFind / MultiLoc-plant

high_taxon_assignment = open("high_taxon_assignment.txt").read().split("\n")
high_taxon_assignment_d = {}
for line in high_taxon_assignment:
	line = line.split("\t")
	try:
		high_taxon_assignment_d[line[0]] = line[1]
	except IndexError:
		print(line, "not found")

outfilepredictions = open(prefix + '-preds.txt','w')
for leaf in taxa:
	query = taxa[leaf]
	if query in preds_d:
		group = query.split("@")[1]
		if group in opisthokonts:
			prediction = ("{}_//_{}".format(preds_d[query]["hectar"], preds_d[query]["ML2ANIMAL"]))
		elif group in otherhetero:
			prediction = ("{}_//_{}".format(preds_d[query]["targetp"], preds_d[query]["ML2ANIMAL"]))
		elif group in primary:
			prediction = ("{}_//_{}".format(preds_d[query]["targetp"], preds_d[query]["ML2PLANT"]))
		elif group in stramenopiles:
			prediction = ("{}_//_{}".format(preds_d[query]["hectar"], preds_d[query]["asafind"]))
		elif group in otherhigher:
			prediction = ("{}_//_{}_//_{}".format(preds_d[query]["targetp"], preds_d[query]["asafind"], preds_d[query]["ML2PLANT"]))
		#zde přidat predikci do taxa dictionary
		prediction = re.sub('[():,]', '', prediction)
		taxa[leaf] = ("'{}__{}'".format(query, prediction))
		outfilepredictions.write("{}\t{}\n".format(query, prediction))
	else:
		genus = query.replace("_", " ")
		genus = genus.split()[0]
		if genus in high_taxon_assignment_d:
			newleaf = ("{}@{}".format(query.split("@")[0], high_taxon_assignment_d[genus]))
			taxa[leaf] = ("'{}'".format(newleaf))


print("prediction summaries written to file: {}-preds.txt".format(prefix))
print("Leaf renaming dictionary ready. Now to tree leaves renaming...")


#taxon renaming in the raxml tree
inTrees = [s for s in os.listdir('.') if s.endswith('.tre') or s.endswith('.treefile')]
for currtree in inTrees:
	if currtree.endswith('-RXM.tre'):
		TREETYPE = "RAxML"
		if "bipartitions" in currtree:
			currtreename = currtree.split(".tre")[0].replace("RAxML_bipartitions.", "")
			tree_line = open(currtree).readline()
			for key in taxa:
				tree_line = tree_line.replace(key, taxa[key])
			with open(currtreename + "-final.tre", "w") as result:
				result.write(tree_line)
		else:
			print("skipping {}, not a tree file".format(currtree))

	elif currtree.endswith('.chain.con.tre'):
		TREETYPE = "PhyloBayes"
		currtreename = currtree.split(".tre")[0].replace("chain.con.", "")
		tree_line = open(currtree).readline()
		for key in taxa:
			tree_line = tree_line.replace(key, taxa[key])
		with open(currtreename + "-PB-final.tre", "w") as result:
			result.write(tree_line)

	elif currtree.endswith('.treefile'):
		TREETYPE = "IQtree"
		currtreename = currtree.split(".tre")[0].replace("trimmed.phy.", "")
		tree_line = open(currtree).readline()
		for key in taxa:
			tree_line = tree_line.replace(key, taxa[key])
		with open(currtreename + "-IQT-final.tre", "w") as result:
			result.write(tree_line)

	elif currtree.endswith('-ASA.tre'):
		TREETYPE = "AsaturA"
		currtreename = currtree.split(".tre")[0]
		tree_line = open(currtree).readline()
		for key in taxa:
			tree_line = tree_line.replace(key, taxa[key])
		with open(currtreename + "-final.tre", "w") as result:
			result.write(tree_line)

	else:
		print(currtree, " is an unrecognized tree file.")

"""
for key in sequencesdictionary.keys():
	code = sequencesdictionary[key].split('_')[0]
	rest = '_'.join(sequencesdictionary[key].split('_')[1:4])
	if key in taxa:
		finalname = '{}_{}'.format(taxa[key], finalpredictions[key])
	elif code in taxa:
		finalname = '{}_{}_{}'.format(taxa[code], rest, finalpredictions[key])
	else:
		finalname = sequencesdictionary[key] + '_' + finalpredictions[key]
	strom_line = strom_line.replace(key, finalname)

with open(prefix + '-fin.tre', 'w') as result:
    result.write(strom_line)
"""
print("Hotovo")
