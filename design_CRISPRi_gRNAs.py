from Bio import Seq
from Bio import SeqIO
import csv
import itertools
import subprocess
import math
from multiprocessing import Pool
import tempfile
import RNA
import argparse
import os
import re
import sys
import random

parser = argparse.ArgumentParser(description='Designs spacers targeting the input genes.')

parser.add_argument('--input_file', '-i', type=str, default="target_genes_sequences.fasta",
					help='Fasta file with sequences of the target genes + promoters. '
					'(Default: target_genes_sequences.fasta)')
					
parser.add_argument('--coordinates', '-c', type=str, default="target_genes_coordinates.txt",
					help='File with coordinates of the target genes. (Default: target_genes_coordinates.txt)')

parser.add_argument('--reference', '-r', type=str, default="genome.fasta", help='Genome sequence file '
					'(Default: genome.fasta).')

parser.add_argument('--spacer_length', '-sl', default="20", type=int,
					help='Length of each spacer. No more than 26 if the -a flag is used. (Default: 20)')
					
parser.add_argument('--promoter_length', '-pl', default="0", type=int,
					help='Length of the promoter region. Keep to 0 if you want to avoid targeting the promoter or '
					'if you\'re not sure where the promoter is. (Default: 0)')
	
parser.add_argument('--non_targeting', '-nt', help='The script designs nontargeting gRNAs if this flag is present. '
					'If only -nt is set, the number of designed nontargeting gRNAs is the largest number between '
					'20 and total number of gRNAs divided by 100. If -nt n is set, the script designs n nontargeting gRNAs.',
					action='store',nargs='?',type=int,default=False,const=9999999)

parser.add_argument('--PAM', '-pam', default="TTV", type=str, help='PAM sequence. (Default: TTV)')

parser.add_argument('--PAM_preference', '-pp', type=str, choices=['template','nontemplate'],
					default="template",	help='Which strand within the coding region should be targeted? (Default: template)')

parser.add_argument('--PAM_orientation', '-po', type=str, choices=["5prime","3prime"],
					default="5prime", help='At which end of the protospacer is the PAM located? (Default: 5prime)')

parser.add_argument('--spacers', '-s', default="4", type=int, help='Number of spacers (including one array if using '
					'the -a flag) to be designed per targeted gene. (Default: 4)')

parser.add_argument('--arrays', '-a', help='The script attempts to design one array per target gene if this flag is present. '
					'This flag can be used only if --PAM is TTV.', action='store_true')

parser.add_argument('--folding', '-f', default=0.2, type=float, help='Minimal accepted correct folding probability of the arrays. '
						 'Scale: 0-1. (Default: 0.2)')					

parser.add_argument('--left_overlap', '-lo', default="atctttgcagtaatttctactgttgtagat", type=str,
					help='Left overhang for cloning grna spacers. (Default: atctttgcagtaatttctactgttgtagat)')	

parser.add_argument('--right_overlap', '-ro', default="ccggcttatcggtcagtttcacctgattta", type=str,
					help='Right overhang for cloning grna spacers. (Default: ccggcttatcggtcagtttcacctgattta)')

parser.add_argument('--processors', '-p', default=1, type=int,
					help='Number of processors. A real impact of using >1 processors is seen only '
					'when also using the -a option. (Default: 1)')

args = parser.parse_args()
processors = args.processors
inname = args.input_file
coord_name = args.coordinates
fold_threshold = args.folding
left_overh = args.left_overlap
right_overh = args.right_overlap
arrays_flag = args.arrays
non_targeting_flag = args.non_targeting
PAM = args.PAM
spacer_length = args.spacer_length
promoter_length = args.promoter_length
promoter_length = args.promoter_length
output_no = args.spacers
reference = args.reference
PAM_preference = args.PAM_preference
PAM_orientation = args.PAM_orientation

if PAM != "TTV" and arrays_flag:
	print('ERROR: you have asked to design arrays (-a) with a PAM different from TTV. Array design is only supported for this PAM.')
	sys.exit()
if PAM != "TTV" and (PAM_preference == "template" or PAM_orientation == "5prime"):
	input('WARNING: you have chosen a PAM different from the default (TTV), but either one or templateh of --PAM_orientation and '
	'--PAM_preference have the default values.\nIf you still want to proceed, press enter.')

if arrays_flag:
	diff_crates = 26 - spacer_length #factor to extend target gene coordinates, 
									#since spacers in arrays to be cloned with the
									#CRATES method (PMID:31270316) need to be 26nt long
	repeat = "GUCUAAGAACUUUAAAUAAUUUCUACUGUUGUAGAU" #sequence of the Cas12a CRISPR repeat
else:
	diff_crates = 0

overlap_window = 10
#window 3' and 5' of spacers in which spacers that fall are defined as overlapping
#NB: the windows sum up! So, if this value is 10, two spacers have to be
#>10+10=20 nt separated to be non-overlapping

stretches = [base*5 for base in ["A","C","G","T"]] #list of base stretches

#create genome_seqs dictionary with sequences of chromosome and plasmids:
genome = SeqIO.parse(reference, 'fasta')
genome_seqs = {}
for record in genome:
	record_id = re.findall('[^ ]+', record.description)[0]
	genome_seqs[record_id] = str(record.seq)

ambiguous_alph = {"N":["A","C","G","T"], "V":["A","C","G"], "H":["A","C","T"], "D":["A","G","T"],
	"B":["C","G","T"], "M":["A","C"], "K":["G","T"], "W":["A","T"], "S":["C","G"], "Y":["C","T"],
	"R":["A","G"], "A": ["A"], "C": ["C"], "G": ["G"], "T": ["T"]}

def remove_ambiguous(sequence):
# returns a list of all possible sequences
# given the ambiguous sequence input (e.g. with N,V, etc. bases)
	return ["".join(i) for i in itertools.product(*[ ambiguous_alph[j] for j in sequence ]) ]

def revcomp(sequence):
# returns reverse complement of sequence
	complement = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N','Y':'R',
	'R':'Y', 'W':'W', 'S':'S', 'K':'M', 'M':'K', 'D':'H', 'H':'D',
	'V':'B', 'B':'V'}
	bases = list(sequence)
	bases.reverse()
	revcomp = ''
	for base in bases:
		if base.islower():
			base = base.upper()
			revcomp += complement[base].lower()
		else:
			revcomp += complement[base]
	return revcomp

def GC_count(string):
# returns %GC of string
	bases = list(string)
	G = bases.count("G")
	C = bases.count("C")
	GC_content = 100*(G+C)/len(string)
	return GC_content

def find_all(string, substr):
# returns a list of the indexes of all occurrencies of substr in string
	listindex=[]
	start = 0
	i = string.find(substr, start)
	while i >= 0:
		listindex.append(i)
		i = string.find(substr, i + 1)
	return listindex

def nontemplate_spacers(PAM_sequence):
# returns a dictionary of start, strand, start_overlap, end_overlap and sequence
# of a spacer targeting the nontemplate (top, +, sense... call it as you wish) 
# strand of the gene_sequence string.
	PAM_indexes = find_all(gene_sequence, revcomp(PAM_sequence))
	spacers_nontemplate = []
	if len(PAM_indexes) != 0:
		for PAM_index in PAM_indexes:
			if PAM_orientation == "5prime":
				start = PAM_index - spacer_length
				end = PAM_index
			else:
				start = PAM_index + PAM_length
				end = PAM_index + PAM_length + spacer_length
			if end > gene_length or start < 0:
				continue
			skip_spacer = filter_spacer(start, end, "nontemplate")
			if skip_spacer:
				continue
			start_overlap = start - overlap_window
			end_overlap = end + overlap_window
			if start_overlap < 0:
				start_overlap = 0 #don't allow indexes to go beyond the target's extremities
			if end_overlap > gene_length:
				end_overlap = gene_length #same as comment above
			start_array = start - diff_crates
			sequence = revcomp(gene_sequence[start_array:end])
			spacers_nontemplate.append({"start":start, "strand":"nontemplate",
			"start_overlap":start_overlap, "end_overlap":end_overlap,"seq":sequence})
	return spacers_nontemplate

def template_spacers(PAM_sequence):
# returns a dictionary of start, strand, start_overlap, end_overlap and sequence
# of a spacer targeting the template (bottom, -, antisense... call it as you wish) 
# strand of the gene_sequence string.
	PAM_indexes = find_all(gene_sequence, PAM_sequence)
	spacers_templatetom = []
	if len(PAM_indexes) != 0:
		for PAM_index in PAM_indexes:
			if PAM_orientation == "5prime":
				start = PAM_index + PAM_length
				end = PAM_index + PAM_length + spacer_length
			else:
				start = PAM_index - spacer_length
				end = PAM_index
			if end > gene_length or start < 0:
				continue
			skip_spacer = filter_spacer(start, end, "template")
			if skip_spacer:
				continue
			start_overlap = start - overlap_window
			end_overlap = end + overlap_window
			if start_overlap < 0:
				start_overlap = 0 #don't allow indexes to go beyond the target's extremities
			if end_overlap > gene_length:
				end_overlap = gene_length #same as comment above
			end_array = end + diff_crates
			sequence = gene_sequence[start:end_array]
			spacers_templatetom.append({"start":start, "strand":"template",
			"start_overlap":start_overlap, "end_overlap":end_overlap,"seq":sequence})
	return spacers_templatetom

def filter_spacer(start, end, strand):
# Checks if the spacer falls within allowed portions of the target gene,
# specified by the --PAM_preference flag. True if it doesn't, False otherwise.
	midpoint = (start+end)/2
	if midpoint in range(0,promoter_length):
		#the spacer targets within the promoter
		return False
	else:
		#the spacer targets downstream of the promoter
		if strand == PAM_preference:
			return False
		else:
			return True

def bestspacers(spacers,n,offt):
# picks the n spacers from the "spacers" list that have the highest scores.
# Scoring system similar to the arrays one: scores the spacer ased on 
# GC content, position within the targeted gene and presence/absence of 
# same-base stretches. The scores of each feature are multiplied to give 
# the score of the whole spacer.
	if offt:#if we are picking spacers from the ones with offtargets, pick from
			#spacers that have the highest number of mismatches
		if len(spacers) <= n:
			return [spacer[0] for spacer in spacers.values()]
		new_spacers_dic = {key:value for key,value in spacers.items() if value[1] == max(value[1] for value in spacers.values())}
		for spacer in new_spacers_dic:
			spacers.pop(spacer)
		new_spacers = [value[0] for key,value in new_spacers_dic.items()]
		best_spacers = bestspacers(new_spacers,n,False)[:]
		new_n = n - len(best_spacers)
		if new_n != 0 and len(spacers) > 0:
			best_spacers.extend(bestspacers(spacers,new_n,True))
		return best_spacers
	best_spacers = []
	if len(spacers) <= n:
		return spacers
	for spacer in spacers:
		#score spacer:
		all_stretches = 0
		seq = spacer["seq"][:spacer_length]
		for x in stretches:
			stretches_no = 0
			indexes = find_all(seq, x)
			stretches_no += len(indexes)
			# count all stretches that are part of a larger stretch as 1
			# (i.e. AAAAAA (six A) counts as only 1, instead of 2)
			if len(indexes) > 1:
				for i in range(len(indexes)-1):
					if indexes[i+1] - indexes[i] == 1:
						stretches_no-= 1
			all_stretches += stretches_no
		GC = GC_count(seq)
		position = spacer["start"]
		if PAM == "TTV":
			#saves the nucleotide 5' of the TTV pam to PAMletter
			try:
				if spacer["strand"] == "template":
					PAMletter = gene_sequence[spacer["start"]-PAM_length-1]
				else:
					PAMletter = revcomp(gene_sequence[spacer["start"] + spacer_length + PAM_length])
			except IndexError:
				#captures when PAMletter would fall outside the gene limits. Sets PAMletter to not T.
				PAMletter = None
			score = calc_spacerScore(GC,all_stretches,position,PAMletter) #score of the spacer
		else:
			score = calc_spacerScore(GC,all_stretches,position,None) #score of the spacer
		spacer["score"] = score
	spacers = sorted(spacers, key=lambda x: (-x['score'], x['start'])) #sort by decreasing score. Spacers with an
																	#equal score are sorted by increasing start 
																	#position, so that those closer to the gene
																	#start are prioritized.
	best_spacers.append(spacers[0]) #append spacer with the highest score
	i = 1
	while len(best_spacers) < n:
		best_spacers.append(spacers[i]) #take spacers from second highest score on
		i += 1
	return best_spacers

def calc_spacerScore(GC,stretch,position,PAMletter):
# calculate score of the spacer
	if 0 <= GC < 15:
		GC_score = 0.1
	elif 15 <= GC < 20:
		GC_score = 0.2
	elif 20 <= GC < 30:
		GC_score = 5
	elif 30 <= GC < 35:
		GC_score = 8	
	elif 35 <= GC < 75:
		GC_score = 10
	elif 75 <= GC < 80:
		GC_score = 8
	elif 80 <= GC < 90:
		GC_score = 4
	elif 90 <= GC:
		GC_score = 0.1
	if stretch == 0:
		stretch_score = 10
	elif stretch == 1:
		stretch_score = 5
	elif stretch > 1:
		stretch_score = 0.1
	if 0 <= position < promoter_length + overlap_window: 
		#if inside the promoter (extended by overlap_window)
		position_score = 10
	else:
		position_score = 1
	PAM_score = 5
	if PAMletter == "T":
		PAM_score = 10
	scores = [GC_score,stretch_score,position_score,PAM_score]
	prod = 1
	#calculate the final score, which is the product of all sub-scores:
	for x in scores:
		prod = prod * x
	return prod

def calc_arrayScore(prob):
# calculate score of the array/spacer based on folding probability
	if .2 <= prob < .3:
		prob_score = 6
	elif .3 <= prob < .4:
		prob_score = 7
	elif .4 <= prob < .5:
		prob_score = 8
	elif prob >= .5:
		prob_score = 10
	return prob_score

def make_arrays(spacers):
# designs Cas12a arrays if spacers has at least 3 spacers. Otherwise, designs simple gRNAs.
	selected = []
	if len(spacers) > 0:
		#pick (up to) n spacers for gRNAs
		n=output_no - 1	
		selected = bestspacers(spacers,n,False)[:]
		#remove these n spacers from the list of available spacers:
		for spacer in selected:
			spacers.remove(spacer)
	# make 1 array (of n spacers) from the "spacers" list and, if not possible, pick single spacers
	if len(spacers) > 2: #try to make array
		#make list of all possible combinations (= position doesn't matter) of n spacers
		combinations = list(itertools.combinations(spacers,n))
		combinations_cp = combinations[:]
		#remove combinations that contain overlapping spacers:
		for comb in combinations_cp:
			comb_sorted = sorted(comb, key = lambda x: x["start"]) #sorts spacers by start index
			i = 0
			while i < 2:
				if comb_sorted[i]["end_overlap"] > comb_sorted[i+1]["start_overlap"]: #if the two spacers overlap
					combinations.remove(comb)
					i = 2
				i += 1
		if len(combinations) > 0:
			#for all non-overlapping combinations, make all permutations and evaluate the folding
			pool = Pool(processors)
			pool_results = pool.map(folding_prob, combinations)
			pool.close()
			pool.join()
			permutations = [] #list of permutations (= position does matter) of n spacers
			for valid_perm in pool_results:
				if len(valid_perm) > 0:
					permutations.extend(valid_perm)
			if len(permutations) > 0:
				permutations = score_arrays(permutations)
				permutations = sorted(permutations, key = lambda x: x[-1], reverse=True) #sort by decreasing array score
				array = list(permutations[0][1]) #permutation with the highest score
				selected.append(array)
				for spacer in array:
					spacers.remove(spacer)
	if len(selected) < output_no and len(spacers) > 1:
		n = output_no - len(selected)
		for spacer in bestspacers(spacers,n,False):
			selected.append(spacer)
	return selected

def score_arrays(arrays):
# gives a score to each array in arrays as follows:
# first it calculates the score of the array based on the folding probability,
# then it scores each spacer of the array individually, based on GC content,
# presence/absence of same-base stretches and position within the targeted gene.
# The scores of each feature are multiplied to give the score of the whole spacer.
# The score of the entire array is then obtained by multiplying the scores of
# all spacers that compose it and the folding probability. The array's score is
# then appended to the array list as last element.
	scored_arrays = []
	for array in arrays:
		prob = array[0]
		#initialize scores list with the folding probability score of the whole array:
		scores = [calc_arrayScore(prob)]
		#score spacers individually and append the score to the scores list:
		for spacer in array[1]:
			all_stretches = 0 #number of stretches
			seq = spacer["seq"]
			shorter_seq = seq[:20]
			for stretch in stretches:
				n = 0
				indexes = find_all(shorter_seq, stretch)
				n += len(indexes)
				# count all stretches that are part of a larger stretch as 1
				# (i.e. AAAAAA (six A) counts as only 1, instead of 2)
				if len(indexes) > 1:
					for i in range(len(indexes)-1):
						if indexes[i+1] - indexes[i] == 1:
							n-= 1
				all_stretches += n
			GC = GC_count(shorter_seq)
			position = spacer["start"] #start index of the spacer
			try:
				if spacer["strand"] == "template":
					PAMletter = gene_sequence[spacer["start"]-PAM_length-1]
				else:
					PAMletter = revcomp(gene_sequence[spacer["start"] + spacer_length + PAM_length])
			except IndexError:
				#captures when PAMletter would fall outside the gene limits. Sets PAMletter to not T.
				PAMletter = None
			score = calc_spacerScore(GC,all_stretches,position,PAMletter) #score of the spacer
			scores.append(score)
		array_score = 1
		#generate the array score by multiplying the scores of all spacers
		#that compose it and the folding probability score:
		for x in scores:
			array_score = array_score * x
		array.append(array_score)
		scored_arrays.append(array)
	return scored_arrays				


def offtarget(thisPAM,thisstrand,spacer_id):
# checks if the offtarget (starting at index; with PAM coordinates start:end)
# is followed by a PAM and if it's outside a targeted gene.
# If so, it's considered as a bona-fide off target
# and the function returns 1, otherwise returns 0.
	global multitargeting
	global offtargeting
	multitarget = False
	if any(thisPAM == PAM for PAM in list_of_PAMs):#check if PAM is a valid PAM
		if "|" in thisgene: #if it's a group of genes having the same sequence
			genes = thisgene.split("|")
			for gene in genes:
				[thisstart,thisend] = positions[scaffold][gene]
				if index in range(thisstart,thisend+1):#if this spacer targets one of these genes,
					return							   #it's not an offtarget
		elif thisgene in positions[scaffold]:
			[thisstart,thisend] = positions[scaffold][thisgene]
		if (thisgene not in positions[scaffold]) or (index not in range(thisstart,thisend+1)): 
		#if doesn't target the current gene (= it' either an offtarget or a multitarget)
			for gene,[i,k] in positions[scaffold].items():
				if gene not in thisgene: #in case that thisgene is a group of genes having the same sequence
					if index in range(i,k+1): #checks if index is inside a target gene
						multitarget = True
						break
			#look for the spacer in the spacers dictionary and remove it
			for spacer in all_spacers_dic[thisgene]:
				if spacer["seq"][:spacer_length] == seq and spacer["strand"] == thisstrand:
					all_spacers_dic[thisgene].remove(spacer) #remove the spacer
					if row[5] != "": #if there are mismatches between the offtarget site and the spacer
						mismatches = 1/(len(row[5].split(","))+1)
					else:
						mismatches = 1
					curr_spacer = spacer
					#append the spacer to the "gene" key of the "offtargeting_spacers" dict.
					#inserts the key into the dict if it doesn't exist:
					if not multitarget: #is true offtarget
						if thisgene not in offtargeting_spacers:
							offtargeting_spacers[thisgene] = {}
						if spacer_id not in offtargeting_spacers[thisgene]:
							offtargeting_spacers[thisgene][thisspacer] = [curr_spacer,0]
						offtargeting_spacers[thisgene][thisspacer][1] += mismatches
						offtargeting += 1
					else: # targets another gene
						if thisgene not in multitargeting_spacers:
							multitargeting_spacers[thisgene] = {}
						if spacer_id not in multitargeting_spacers[thisgene]:
							multitargeting_spacers[thisgene][thisspacer] = [curr_spacer,0]
						multitargeting_spacers[thisgene][thisspacer][1] += mismatches
						multitargeting += 1
					break

def folding_prob(combination):
# returns a list of all permutations (= position DOES matter) of 3 spacers 
# whose array has a probability of folding correctly >= fold_threshold
	valid_perms = []
	permutations = list(itertools.permutations(combination,3))
	for perm in permutations:
		spacers_seqs = [x["seq"] for x in perm]
		prob = folding(spacers_seqs) #folding probability of the array
		print(prob)
		if prob >= fold_threshold:
			valid_perms.append([prob,perm,combination])
	return valid_perms

def folding(spacers):
# calculate the probability of the array folding in the desired structure
	array = ""#array sequence
	repeat_c = "...................xx(((((....)))))." #desired folding of the array
													  #with stem-loop at the 3' end
	array_c = "" #array dot-bracket constraint
	for spacer in spacers:
		array = array + repeat + spacer
		spacer_c = "."*len(spacer)
		array_c = array_c + repeat_c + spacer_c
	if "CGTCTC" in array or revcomp("CGTCTC") in array:
		#remove arrays that contain a BsmBI site by giving them
		#a negative folding probability
		return -1
	R=0.00198 #Boltzmann constant
	T=310.15 #37 celsius in Kelvin
	# create a fold_compound object for the current sequence:
	fc = RNA.fold_compound(array)
	# compute the MFE and corresponding structure:
	(mfe_struct, mfe) = fc.mfe()
	# re-scale the Boltzmann factors for partition function computation
	fc.exp_params_rescale(mfe)
	# compute partition function:
	(pp, Fu) = fc.pf()
	# re-create a fold_compound object:
	fc = RNA.fold_compound(array)
	# add hard constraint from dot-bracket notation:
	fc.hc_add_from_db(array_c)
	# compute the new MFE and corresponding structure with the constraints:
	(mfe_struct, mfeNew) = fc.mfe()
	# re-scale the Boltzmann factors for partition function computation:
	fc.exp_params_rescale(mfeNew)
	# compute partition function:
	(pp, Fc) = fc.pf()
	#calculate probability of folding within the constraints:
	array_prob=math.exp((Fu-Fc)/(R*T))
	return array_prob

def random_string(length):
# makes a random sequence of length bases
	return ''.join(random.choice(['A','C','G','T']) for i in range(length))

#create a dictionary of start and end positions of the targeted genes
#of structure scaffold:[[start1,end1],[start2,end2]]:
positions = {}
genes_coord = open(coord_name,mode="r")
for line in genes_coord:
	if not line.rstrip():
		continue
	line = line.split("\t")
	scaffold = line[0]
	gene_id = line[1]
	start = int(line[2])
	end = int(line[3])
	if scaffold not in positions:
		positions[scaffold] = {}
	positions[scaffold][gene_id] = [start,end]

#make some variables used later:
gRNAs_per_gene = {} #for storing gRNAs selected for each gene
arrays_per_gene = {} #for storing arrays selected for each gene
gene_list = SeqIO.parse(inname, 'fasta')
gene_seqs = {} #dic of target_gene_id:target_gene_sequence
for gene in gene_list:
	gene_seqs[gene.id] = str(gene.seq)
all_spacers_dic = {} #all possible spacers for each PAM, before any filtering
all_PAMs = [] #list of all PAMs (without IUPAC ambiguity codes)
all_revcomp_PAMs = []
PAM_list = remove_ambiguous(PAM)
all_PAMs.extend(PAM_list)
revcomp_PAM_list = [revcomp(PAM) for PAM in PAM_list]
all_revcomp_PAMs.extend(revcomp_PAM_list)
PAM_length = len(PAM)
gene_IDs = []

###Find all spacers and filter offtargets:
print("Looking for spacers and offtargets...\n")
##find all spacers:
for gene_id, gene_sequence in gene_seqs.items():
	gene_IDs.append(gene_id)
	gene_length = len(gene_sequence)
	spacers_nontemplate = []
	spacers_templatetom = []
	for PAM_sequence in PAM_list:
		spacers_nontemplate.extend(nontemplate_spacers(PAM_sequence))
		spacers_templatetom.extend(template_spacers(PAM_sequence))
		all_spacers = spacers_nontemplate + spacers_templatetom
	if gene_id not in all_spacers_dic:
		all_spacers_dic[gene_id] = [] 
	all_spacers_dic[gene_id] = all_spacers

##Run bowtie to find offtargeting spacers:
#make input file for bowtie
bowtie_in = tempfile.NamedTemporaryFile(mode='w+') #input file for bowtie
for gene_id, selectedspacers in all_spacers_dic.items():
	i = 0
	for spacer in selectedspacers:
		seq = spacer["seq"][:spacer_length] #look only for offtargets in the spacer_length part of the spacer
		strand = spacer["strand"]
		bowtie_in.write(">")
		bowtie_in.write("{0} {1} {2}".format(gene_id, i,strand)) #ID of the fasta entry: targeted gene + i + targetedstrand
		bowtie_in.write("\n")
		bowtie_in.write(seq) #spacer sequence
		bowtie_in.write("\n")
		i += 1
bowtie_in.seek(0)

#make bowtie indexes if not present:
if not os.path.exists("bowtie_files"):
	os.mkdir("bowtie_files")
if not os.path.isfile('bowtie_files/*.ebwt'):
	cmd = ("bowtie-build --threads {} -r -q -f {} bowtie_files/genome_index".format(processors,reference))
	try:
		subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT,cwd=None)
	except subprocess.CalledProcessError as e:
		print("ERROR: Creating the bowtie index with the command: \n{} \nfailed with this message: \n{}".format(e.cmd,e.output.decode('utf-8')))
		sys.exit()

#run bowtie:
bowtie_out = tempfile.NamedTemporaryFile(mode='w+')
offtargeting = 0 #no. of spacers discarded because of off-targeting
multitargeting = 0
cmd = ("bowtie --threads {0} -v 2 -a bowtie_files/genome_index --suppress 6,7 "
"-f {1} > {2}".format(processors,bowtie_in.name,bowtie_out.name)) #matches with up to 2 mismatches (included)
try:
	subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT,cwd=None)
except subprocess.CalledProcessError as e:
	print("ERROR: Bowtie run with the command: \n{} \nfailed with this message: \n{}".format(e.cmd,e.output.decode('utf-8')))
	sys.exit()

#read bowtie results
offtargeting_spacers = {}
multitargeting_spacers = {}
csv_reader = csv.reader(bowtie_out, delimiter='\t')
rowlist = []
for row in csv_reader:
	rowlist.append(row)
rowlist = sorted(rowlist, key = lambda x: x[0]) # sort alphabetically by spacer sequence, so that
												# we can stop at the first offtarget of a spacer
for row in rowlist:
	thisspacer = row[0]
	(thisgene, _, thisstrand) = thisspacer.split(" ")
	strand = row[1]
	scaffold = row[2] #chromosome or plasmid ID
	index = int(row[3]) #starting index of off-target
	seq = row[4]
	if strand == "-":
		seq = revcomp(seq)
	# get start and end positions of the PAM of the off-target and check if it's a valid PAM
	# if it is, the off-target is a true one
	if PAM_orientation == "5prime":
		if strand == "+":
			start = index - PAM_length
			end = index
			thisPAM = genome_seqs[scaffold][start:end]
		else:
			start = index + spacer_length
			end = index + spacer_length + PAM_length
			thisPAM = revcomp(genome_seqs[scaffold][start:end])
	else:
		if strand == "+":
			start = index + spacer_length
			end = index + spacer_length + PAM_length
			thisPAM = genome_seqs[scaffold][start:end]
		else:
			start = index - PAM_length
			end = index
			thisPAM = revcomp(genome_seqs[scaffold][start:end])
	offtarget(thisPAM,thisstrand,thisspacer)
bowtie_in.close()
bowtie_out.close()

###design gRNA/arrays:
if arrays_flag == True:
	print("Designing gRNAs and arrays (this may take long)...\n")
	print("Currently working on gene:")
else:
	print("Designing gRNAs...")
	
for gene_id, spacers in all_spacers_dic.items():
	gene_sequence = gene_seqs[gene_id]
	gene_length = len(gene_sequence)
	offt = False
	if arrays_flag == True:
		print(gene_id)
		#remove spacers whose spacer_length+diff_crates long sequence
		#falls outside the target gene:
		spacers = [x for x in spacers if len(x["seq"]) == spacer_length+diff_crates]
		selected = make_arrays(spacers)
	else:
		selected = bestspacers(spacers,output_no,False)[:]
	if len(selected) < output_no:
		n = output_no - len(selected)
		#if there are not enough (== output_no) selected arrays/spacers,
		#pick spacers from the ones that have multipletargets
		if gene_id in multitargeting_spacers:
			selected.extend(bestspacers(multitargeting_spacers[gene_id],n,True))
	if len(selected) < output_no:
		n = output_no - len(selected)
		#if there are not enough (== output_no) selected arrays/spacers,
		#pick spacers from the ones that have off-targets
		if gene_id in offtargeting_spacers:
			selected.extend(bestspacers(offtargeting_spacers[gene_id],n,True))
	for j in selected:
		if len(j) == 3:
			#it's an array
			arrays_per_gene[gene_id] = [x["seq"] for x in j]
		else:
			#it's a gRNA
			seq = j["seq"]
			gRNAs_per_gene.setdefault(gene_id,[]).append(seq[:spacer_length])

#finds spacers that target more genes. These spacers survived the offtargeting filter
#because not enough spacers targeting the gene could be designed.
rev_dict = {}
for gene,spacers in gRNAs_per_gene.items():
	for spacer in spacers:
		rev_dict.setdefault(spacer, []).append(gene)

#group these spacers under the name gene1_gene2:
for spacer,genes in rev_dict.items():
	if len(genes) > 1:
		for gene in genes:
			gRNAs_per_gene[gene].remove(spacer)
		gRNAs_per_gene.setdefault("_".join(genes), []).append(spacer)

#remove genes that have all spacers in common with another:
gRNAs_per_gene = {key:val for key, val in gRNAs_per_gene.items() if val != []}
new_gene_IDs = gene_IDs + [x for x in gRNAs_per_gene if x not in gene_IDs]

gRNA_counter = sum([len(x) for x in gRNAs_per_gene.values()]) #counts how many gRNAs were designed

if non_targeting_flag:
	##make non-targeting gRNAs:
	#find out how many to make:
	if non_targeting_flag == 9999999: #default value of -nt
		min_output_size = max(20,int(gRNA_counter/100))
	else:
		min_output_size = non_targeting_flag
	print("\nDesigning {} non targeting gRNAs...".format(min_output_size))
	##make the gRNAs:
	nt_spacers = []
	while len(nt_spacers) < min_output_size:
		with tempfile.NamedTemporaryFile(mode='w+') as outf:
			strings = []
			#make 1000 random spacers and filter by %GC:
			for n in range(1000):
				string = random_string(spacer_length)
				if 30<= GC_count(string) <= 70:
					strings.append(string)
			#remove targeting ones
			with tempfile.NamedTemporaryFile(mode='w+') as tmp_in:
				i = 0
				for gRNA in strings:
					tmp_in.write(">nt" + str(i) + "\n")
					tmp_in.write(gRNA + "\n")
					i += 1
				tmp_in.flush()
				#run bowtie:
				cmd = ("bowtie --threads {0} -v 3 -a bowtie_files/genome_index --suppress 6,7 "
				"-f {1} --un {2}".format(processors,tmp_in.name,outf.name)) #matches with up to 3 mismatches
				subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
			outf.flush()
			out_parse = SeqIO.parse(outf, 'fasta')
			for record in out_parse:
				nt_spacers.append(str(record.seq))
			nt_spacers = list(set(nt_spacers)) # removes duplicates
	tmp_in.close()

subprocess.call("rm -r bowtie_files",shell=True) #remove bowtie folder

references_out = open("references.fasta",mode="w")

#Design oligonucleotides:
if arrays_flag and len(arrays_per_gene) > 0:
	#make and save array oligonucleotides:
	junctions = ["CTAA", "GACA", "GCAC", "AATC", "GTAA", "TGAA", "ATTA", "CCAG",
	"AGGA", "ACAA", "TAGA", "CGGA", "CATA", "CAGC", "AAGT", "CTCC", "AGAT", "ACCA",
	"AGTG", "GGTA", "GCGA", "AAAA", "ATGA"]
	start_junction = "CCTC"
	end_junction = "AACG"
	repeat_dna = "gtctaagaactttaaataatttctactgttgtagat"
	rev_repeat_dna = revcomp(repeat_dna)
	array_pools = int(len(arrays_per_gene)/len(junctions))+1 #number of separate pools of array oligonucleotides
	#make pools of array oligonucleotides:
	n = 0
	for pool_number in range(1,array_pools+1):
		i = 0
		out_arrays = open("array_oligo_pool_{0}.txt".format(pool_number),mode="w")
		while i < len(junctions) and n < len(arrays_per_gene):
			gene = list(arrays_per_gene)[n]#gene name
			array = arrays_per_gene[gene]
			junction = junctions[i]
			revc_junction = revcomp(junction)
			fw1 = start_junction + repeat_dna + array[0] + junction #forward 1 oligo
			rev1 = revcomp(array[0]) + rev_repeat_dna #reverse 1 oligo
			fw2 = repeat_dna + array[1] #forward 2 oligo
			rev2 = revc_junction + revcomp(array[1]) + rev_repeat_dna + revc_junction #etc.
			fw3 = junction + repeat_dna + array[2]
			rev3 = revcomp(end_junction) + revcomp(array[2]) + rev_repeat_dna
			n += 1
			out_arrays.write("\n".join([fw1,rev1,fw2,rev2,fw3,rev3]))
			out_arrays.write("\n")
			construct = "".join([fw1,fw2,fw3,end_junction])
			i += 1
			#replace the value with the actual array sequence:
			arrays_per_gene[gene] = construct

#save the gRNA oligonucleotides:
out_gRNAs = open("gRNA_oligo_pool.txt",mode="w") #file with gRNA oligonucleotides
for gene in new_gene_IDs:
	i = 1
	if gene in gRNAs_per_gene:
		for spacer in gRNAs_per_gene[gene]:
			oligo = left_overh + spacer + right_overh
			out_gRNAs.write(oligo)
			out_gRNAs.write("\n")
			line = ">" + gene + "_grna_" + str(i) + "\n"
			references_out.write(line)
			line = oligo + "\n"
			references_out.write(line)
			i += 1
	if arrays_flag:
		if gene in arrays_per_gene:
			line = ">" + gene + "_array" + "\n"
			references_out.write(line)
			line = arrays_per_gene[gene] + "\n"
			references_out.write(line)

if non_targeting_flag:	
	i = 1
	for nt_gRNA in nt_spacers[:min_output_size]:
		oligo = left_overh + nt_gRNA[:spacer_length] + right_overh
		out_gRNAs.write(oligo)
		out_gRNAs.write("\n")
		line = "non_targeting_" + str(i) + "\n"
		references_out.write(line)
		line = oligo + "\n"
		references_out.write(line)
		i += 1
out_gRNAs.close()


print("\nDone!\n")

if non_targeting_flag:
	msg = "Total number of designed gRNAs: {0} (+ {1} nontargeting).\n".format(gRNA_counter,min_output_size)
else:
	msg = "Total number of designed gRNAs: {0}.\n".format(gRNA_counter)
print(msg)

nogrna_names = [x for x in new_gene_IDs if x not in gRNAs_per_gene]
nogRNA_counter = len(nogrna_names)

if nogRNA_counter > 0:
	print("For {0} gene(s) no gRNAs could be designed: {1}.\n".format(nogRNA_counter,", ".join(nogrna_names)))
else:
	print("At least one gRNA per gene could be designed.\n")
	
if arrays_flag == True:
	array_genes = len(arrays_per_gene)
	if array_genes == 0: #if the list is empty, therefore no arrays have been designed
		print("No arrays could be designed for the input targets.")
	else:
		print("\nNumber of array-cloning reactions:", array_pools)
		print("\nTarget genes with:\n0 arrays: {0}\n1 array: {1}".format(len(gene_IDs)-array_genes, array_genes))