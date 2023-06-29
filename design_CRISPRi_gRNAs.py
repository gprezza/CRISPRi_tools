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
from datetime import datetime
import sys
import random

startTime = datetime.now()

parser = argparse.ArgumentParser(description='Designs spacers '
								'targeting the input genes.')

parser.add_argument('--processors', '-p', default=1, type=int,
					help='Number of processors (Default: 1).')

parser.add_argument('--in_name', '-i', type=str, default="target_genes",
					help='Input base name of the fasta file with sequences of the target genes.'
					'(Default: target_genes.)')
					
parser.add_argument('--coordinates', '-c', type=str,
					help='File with coordinates of the target genes. (Default: "--in_name + _coordinates.txt".)')

parser.add_argument('--folding', '-f', default=0.2, type=float,
					help='Minimal accepted correct folding probability of the arrays. '
						 'Scale: 0-1. (Default: 0.2).')					

parser.add_argument('--left', '-lo', default="atctttgcagtaatttctactgttgtagat", type=str,
					help='Left overhang for cloning mature spacers. (Default: atctttgcagtaatttctactgttgtagat).')	

parser.add_argument('--right', '-ro', default="ccggcttatcggtcagtttcacctgattta", type=str,
					help='Right overhang for cloning mature spacers. (Default: ccggcttatcggtcagtttcacctgattta).')

parser.add_argument('--arrays', '-a', help='The script attempts to design one array per target gene if this flag is present.',
	action='store_true')

parser.add_argument('--non_targeting', '-nt', help='The script designs nontargeting gRNAs if this flag is present.',
	action='store_true')

parser.add_argument('--PAM', '-pam', default="TTV", type=str,
					help='PAM sequence. (Default: TTV).')

parser.add_argument('--spacer_length', '-sl', default="20", type=int,
					help='Length of each spacer. No more than 26 if the -a flag is used. (Default: 20).')
					
parser.add_argument('--promoter_length', '-pl', default="50", type=int,
					help='Length of the promoter region. (Default: 50).')
					
parser.add_argument('--mutants', '-m', default="4", type=int,
					help='Number of spacers (including one array if using the -a flag) to be designed per targeted gene. (Default: 4).')

parser.add_argument('--reference', '-r', type=str, required=True,
					help='Folder containing the genome sequence(s).')


args = parser.parse_args()
processors = args.processors
inname = args.in_name
if args.coordinates is None:
    coord_name = inname
else:
	coord_name = args.coordinates
fold_threshold = args.folding
left_overh = args.left
right_overh = args.right
arrays_flag = args.arrays
non_targeting_flag = args.non_targeting
PAM = args.PAM
spacer_length = args.spacer_length
promoter_length = args.promoter_length
promoter_length = args.promoter_length
output_no = args.mutants
reference = args.reference

#Parameters:
diff_crates = 26 - spacer_length #factor to extend target gene coordinates, 
								 #since spacers in arrays to be cloned with the
							     #CRATES method (PMID:31270316) need to be 26nt long

overlap_window = 10
#window 3' and 5' of spacers in which spacers that fall are defined as overlapping
#NB: the windows sum up! So, if this value is 10, two spacers have to be
#>10+10=20 nt separated to be non-overlapping

repeat = "GUCUAAGAACUUUAAAUAAUUUCUACUGUUGUAGAU" #sequence of the Cas12a CRISPR repeat

stretches = [base*5 for base in ["A","C","G","T"]] #list of base stretches

#make a unique multifasta file with all sequences present in the folder specified by the -r flag:
genome_fa_file = open("{0}/genome_sequence.fa".format(reference),mode= "w+")
for filename in os.listdir(reference):
	if filename != "genome_sequence.fa":
		_file = open("{0}/{1}".format(reference,filename),mode= "r+")
		for line in _file:
			genome_fa_file.write(line)
		genome_fa_file.write("\n")
genome_fa_file.seek(0)

#create genome_seqs dictionary with sequences of chromosome and plasmids:
genome = SeqIO.parse(genome_fa_file, 'fasta')
genome_seqs = {}
for record in genome:
	record_id = re.findall('[^ ]+', record.description)[0]
	genome_seqs[record_id] = str(record.seq)

ambiguous_alph = {"N":["A","C","G","T"], "V":["A","C","G"], "H":["A","C","T"], "D":["A","G","T"],
	"B":["C","G","T"], "M":["A","C"], "K":["G","T"], "W":["A","T"], "S":["C","G"], "Y":["C","T"],
	"R":["A","G"], "A": ["A"], "C": ["C"], "G": ["G"], "T": ["T"]}

def remove_ambiguous(sequence):
	#returns a list of all possible sequences
	#given the ambiguous sequence input (e.g. with N,V, etc. bases)
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

def top_spacers(search_area, PAM_sequence):
# returns a list of lists [a,b,c,d,e,f] where a and b are the indexes
# of a spacer targeting the non template (top, +, sense... call it as you wish) 
# strand of the search_area string; c is top (to indicate the strand), 
# d and e are the spacer indexes + overlap_window and f is the sequence of the spacer.
	PAM_indexes = find_all(search_area, revcomp(PAM_sequence))
	spacers_top = []
	if len(PAM_indexes) != 0:
		for PAM_index in PAM_indexes:
			#check that the index and the subsequent spacer fall within the gene limits
			if PAM_index >= spacer_length and PAM_index <= gene_length - diff_crates:
				start = PAM_index - spacer_length
				if start < diff_crates:
					continue
				end = PAM_index
				start_overlap = start - overlap_window
				end_overlap = end + overlap_window
				if start_overlap < 0:
					start_overlap = 0 #don't allow indexes to go beyond the target's extremities
				if end_overlap > gene_length:
					end_overlap = gene_length #same as comment above
				start_array = start - diff_crates
				sequence = revcomp(search_area[start_array:end])
				PAMletter = revcomp(gene_sequence[PAM_index+PAM_length])
				spacers_top.append([start, PAMletter, "top", start_overlap, end_overlap,sequence])
	return spacers_top

def bottom_spacers(search_area, PAM_sequence):
# returns a list of lists [a,b,c,d,e,f] where a and b are the indexes
# of a spacer targeting the template (bottom, -, antisense... call it as you wish) 
# strand of the search_area string; c is 1 (to indicate the strand),
# d and e are the spacer indexes + overlap_window and f is the sequence of the spacer.
	PAM_indexes = find_all(search_area, PAM_sequence)
	spacers_bottom = []
	if len(PAM_indexes) != 0:
		for PAM_index in PAM_indexes:
			#check that the index and the subsequent spacer fall within the gene limits
			if PAM_index <= gene_length - spacer_length - PAM_length and PAM_index >= diff_crates:
				end = PAM_index + PAM_length + spacer_length
				if end > gene_length - diff_crates:
					continue
				start = PAM_index + PAM_length
				start_overlap = start - overlap_window
				end_overlap = end + overlap_window
				if start_overlap < 0:
					start_overlap = 0 #don't allow indexes to go beyond the target's extremities
				if end_overlap > gene_length:
					end_overlap = gene_length #same as comment above
				end_array = end + diff_crates
				sequence = search_area[start:end_array]
				PAMletter = search_area[PAM_index-1]
				spacers_bottom.append([start, PAMletter, "bot", start_overlap, end_overlap,sequence])
	return spacers_bottom

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
			seq = spacer[5]
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
			position = spacer[0] #start index of the spacer
			PAMletter = spacer[1]
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

def bestspacers(spacers,n,offt):
# picks the n spacers from the "spacers" list that have the highest scores.
# Scoring system similar to the arrays one: scores the spacer ased on 
# GC content, position within the targeted gene and presence/absence of 
# same-base stretches. The scores of each feature are multiplied to give 
# the score of the whole spacer.
	if offt:#if we are picking spacers from the ones with offtargets, just pick from
			#spacers that have the highest number of mismatches
		#sort by decreasing number of mismatches with the off-targeting sites:
		#ordered_spacers = [k for k, v in sorted(spacers.items(), key = lambda item : item[1][1])]
		#print("ppp",[key for key,value in spacers.items() if value[1] == max(value[1] for value in spacers.values())])
		if len(spacers) <= n:
			return [spacer[0] for spacer in spacers.values()]
		new_spacers_dic = {key:value for key,value in spacers.items() if value[1] == max(value[1] for value in spacers.values())}
		for spacer in new_spacers_dic:
			spacers.pop(spacer)
		new_spacers = [value[0] for key,value in new_spacers_dic.items()]
		best_spacers = bestspacers(new_spacers,n,False)
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
		seq = spacer[5][:20]
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
		position = spacer[0]
		PAMletter = spacer[1]
		score = calc_spacerScore(GC,all_stretches,position,PAMletter) #score of the spacer
		spacer.append(score)
	spacers = sorted(spacers, key = lambda x: x[-1],reverse=True) #sort by decreasing score
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
	if 0 <= position < promoter_length + overlap_window + diff_crates: 
		#if inside the promoter (extended by overlap_window)
		position_score = 10
	else:
		#calculate position as percentage inside the coding sequence
		n = promoter_length + overlap_window + diff_crates
		position = (position-n)*100/(gene_length-n) 
		if 0 <= position <= 50:
			position_score = 9
		elif 50 < position <= 60:
			position_score = 7
		elif 60 < position <= 70:
			position_score = 4
		elif 70 < position <= 80:
			position_score = 3
		elif 80 < position <= 90:
			position_score = 2
		elif 90 < position <= 100:
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
# calculate score of the array/spacer
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
	selected = []
	if len(spacers) > 0:
		#pick (up to) n spacers for gRNAs
		n=output_no - 1	
		selected = bestspacers(spacers,n,False)
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
			comb_sorted = sorted(comb, key = lambda x: x[0]) #sorts spacers by start index
			i = 0
			while i < 2:
				if comb_sorted[i][4] > comb_sorted[i+1][3]: #if the two spacers overlap
					combinations.remove(comb)
					i = 2
				i += 1
		if len(combinations) > 0:
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

def offtarget(list_of_PAMs,target_strand,spacer_id):
#checks if the offtarget (starting at index; with PAM coordinates start:end)
#is followed by a PAM and if it's outside a targeted gene.
#If so, it's considered as a bona-fide off target
#and the function returns 1, otherwise returns 0.
	global multitargeting
	global offtargeting
	n = 0
	multitarget = False
	if any(genome_seqs[scaffold][start:end] == PAM for PAM in list_of_PAMs):#check if PAM is a valid PAM
		if "|" in thisgene: #if it's a group of genes having the same sequence
			genes = thisgene.split("|")
			for gene in genes:
				[thisstart,thisend] = positions[scaffold][gene]
				if index in range(thisstart,thisend+1):#if this spacer targets one of these genes,
					return							   #it's not an offtarget
		else:
			[thisstart,thisend] = positions[scaffold][thisgene]
		if not index in range(thisstart,thisend+1): #if doesn't target the current gene (= it' either an offtarget
													#or a multitarget)
			for gene,[i,k] in positions[scaffold].items():
				if gene not in thisgene: #in case that thisgene is a group of genes having the same sequence
					if index in range(i,k+1): #checks if index is inside a target gene
						multitarget = True
						break
			#look for the spacer in the spacers dictionary and remove it
			for spacer in all_spacers_dic[thisgene]:
				if spacer[5][:20] == seq and spacer[2] == target_strand:
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
	#returns a list of all permutations (= position DOES matter) of 3 spacers 
	#whose array has a probability of folding correctly >= fold_threshold
	valid_perms = []
	permutations = list(itertools.permutations(combination,3))
	for perm in permutations:
		spacers_seqs = [x[5] for x in perm]
		prob = folding(spacers_seqs) #folding probability of the array
		if prob >= fold_threshold:
			valid_perms.append([prob,perm,combination])
	return valid_perms

def folding(spacers):
	#calculate the probability of the array folding in the desired structure
	array = ""#array sequence
	repeat_c = "...................xx(((((....)))))." #desired folding of the array
													  #with stem-loop at the 3' end
	array_c = "" #array dot-bracket constraint
	for spacer in spacers:
		array = array + repeat + spacer
		spacer_c = "."*len(spacer)
		array_c = array_c + repeat_c + spacer_c
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
	return ''.join(random.choice(['A','C','G','T']) for i in range(length))
		
#------------- main:

#create a dictionary of start and end positions of the targeted genes
#of structure scaffold:[[start1,end1],[start2,end2]]:
positions = {}
if args.coordinates is None:
	genes_coord = open("{}_coordinates.txt".format(coord_name),mode="r")
else:
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

#create a dictionary of start and end positions of the targeted genes
#of structure scaffold:[[start1,end1],[start2,end2]]:
selected_per_gene = {} #for storing selected spacers/arrays for each PAM
gene_list = SeqIO.parse("{}_sequences.fasta".format(inname), 'fasta')
gene_seqs = {} #dic of target_gene_id:target_gene_sequence
for gene in gene_list:
	gene_seqs[gene.id] = str(gene.seq)
	selected_per_gene[gene.id] = []
all_spacers_dic = {} #all possible spacers for each PAM, before any filtering
all_PAMs = [] #list of all PAMs (without IUPAC ambiguity codes)
all_revcomp_PAMs = []
PAM_list = remove_ambiguous(PAM)
all_PAMs.extend(PAM_list)
revcomp_PAM_list = [revcomp(PAM) for PAM in PAM_list]
all_revcomp_PAMs.extend(revcomp_PAM_list)
PAM_length = len(PAM)
#index until which to extend the promoter, in order to include the spacers
#that target up to spacer_length/2 nt after the TSS:
lo = promoter_length + int(spacer_length/2) + PAM_length + diff_crates
for gene_id, gene_sequence in gene_seqs.items():
	#print(gene_id)
	gene_length = len(gene_sequence)
	longer_promoter = gene_sequence[:lo] #to include spacers targeting up to spacer_length/2 nt after TSS
	spacers_top = []
	spacers_bottom = []
	for PAM_sequence in PAM_list:
		#promoter:
		spacers_top.extend(top_spacers(longer_promoter, PAM_sequence))
		#rest of gene:
		spacers_bottom.extend(bottom_spacers(gene_sequence, PAM_sequence))
		all_spacers = spacers_top + spacers_bottom
	if gene_id not in all_spacers_dic:
		all_spacers_dic[gene_id] = [] 
	all_spacers_dic[gene_id] = all_spacers

previous_gene = ""
bowtie_in = tempfile.NamedTemporaryFile(mode='w+') #input file for bowtie
previous_gene = ""
for gene_id, selectedspacers in all_spacers_dic.items():
	i = 0
	for spacer in selectedspacers:
		seq = spacer[5][:spacer_length] #look only for offtargets in the spacer_length part of the spacer
		strand = spacer[2]
		bowtie_in.write(">")
		bowtie_in.write("{0}_{1}_{2}".format(gene_id, i,strand)) #ID of the fasta entry: targeted gene + _i + _targetedstrand
		bowtie_in.write("\n")
		bowtie_in.write(seq) #spacer sequence
		bowtie_in.write("\n")
		i += 1
bowtie_in.seek(0)

#make bowtie indexes if not present:
if not os.path.exists("bowtie_files"):
	os.mkdir("bowtie_files")
if not os.path.isfile('bowtie_files/*.ebwt'):
	cmd = ("bowtie-build --threads {} -r -q -f {}/genome_sequence.fa bowtie_files/genome_index".format(processors,reference))
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
"-f {1} > {2}".format(processors,bowtie_in.name,bowtie_out.name)) #matches with up to 3 mismatches
try:
	subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT,cwd=None)
except subprocess.CalledProcessError as e:
	print("ERROR: Bowtie run with the command: \n{} \nfailed with this message: \n{}".format(e.cmd,e.output.decode('utf-8')))
	sys.exit()

offtargeting_spacers = {}
multitargeting_spacers = {}

csv_reader = csv.reader(bowtie_out, delimiter='\t')
rowlist = []
for row in csv_reader:
	rowlist.append(row)
rowlist = sorted(rowlist, key = lambda x: x[0]) # sort alphabetically by spacer sequence, so that
												# we can stop at the first offtarget of a spacer
for row in rowlist:
	#print(row)
	thisspacer = row[0]
	(thisgene, _, target_strand) = thisspacer.split("_")
	strand = row[1]
	scaffold = row[2] #chromosome or plasmid ID
	index = int(row[3]) #starting index of off-target
	if strand == "-":
		seq = revcomp(row[4])
		start = index + spacer_length
		end = index + spacer_length + PAM_length
		offtarget(all_revcomp_PAMs,target_strand,thisspacer)
	else:
		seq = row[4]
		start = index - PAM_length
		end = index
		offtarget(all_PAMs,target_strand,thisspacer)
bowtie_in.close()
bowtie_out.close()

total_mature = 0
for gene_id, spacers in all_spacers_dic.items():
	print("\n", gene_id,"\n")
	gene_sequence = gene_seqs[gene_id][diff_crates:-diff_crates]
	gene_length = len(gene_sequence)
	offt = False
	if arrays_flag == True:
		selected = make_arrays(spacers)
	else:
		selected = bestspacers(spacers,output_no,False)
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
	out = []
	for j in selected:
		if len(j) == 3:
			out.append([x[5] for x in j])
		else:
			seq = j[5]
			out.append(seq[:spacer_length])
	selected_per_gene[gene_id] = out
	total_mature += len(out)
	print(out)

if non_targeting_flag:
	#make non-targeting gRNAs:

	print("making nt grnas")
	PAM_no = len(all_PAMs)

	generated_guides = 50000
	min_output_size = max(20,int(total_mature/100))
	limit = 10 #base position up to which the mismathces are checked

	b = 1
	nt_spacers = []

	with tempfile.NamedTemporaryFile(mode='w+') as outf:
		while len(nt_spacers) < min_output_size:
			strings = []
			for i in range(generated_guides):
				string = random_string(spacer_length)
				if 30<= GC_count(string) <= 70:
					strings.append(string)
			with tempfile.NamedTemporaryFile(mode='w+') as tmp_in:
				i = 0
				for gRNA in strings:
					tmp_in.write(">" + str(i) + '_' + "lol" + "\n")
					tmp_in.write(gRNA + "\n")
					i += 1
				tmp_in.flush()
				#run bowtie:
				cmd = ("bowtie --threads {0} -v 1 -3 8 -a bowtie_files/genome_index --suppress 6,7 "
				"-f {1} --un {2}".format(processors,tmp_in.name,outf.name)) #matches with up to 3 mismatches
				subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
			outf.flush()
			out_parse = SeqIO.parse(outf, 'fasta')
			nontargeting = [[0, ''] for j in range(generated_guides)]
			for record in out_parse:
				nt_spacers.append(str(record.seq))
			nt_spacers = list(set(nt_spacers)) # removes duplicates
			#print("Cycle run ", b, "times... I found", len(nt_spacers), "non-targeting gRNAs:")
			b += 1
	#print("Non targeting gRNAs:",nt_spacers[:min_output_size])

#make oligos:

arrays = []
mature = {}
noarray = 0
onearray = 0
for gene, gRNAs in selected_per_gene.items():
	#print(gene)
	#print(gRNAs)
	n = 0
	for group in gRNAs:
		if len(group) == 3:
			arrays.append([gene,group])
			n=1
		else:
			if gene in mature:
				mature[gene].append(group)
			else:
				mature[gene] = [group]
	if n == 0:
		noarray += 1
	else:
		onearray += 1

#finds spacers that target more sRNAs:
rev_dict = {}
for gene,spacers in mature.items():
	for spacer in spacers:
		rev_dict.setdefault(spacer, set()).add(gene)

#group these spacers under the name srna1_srna2:
for spacer,genes in rev_dict.items():
	if len(genes) > 1:
		for gene in genes:
			mature[gene].remove(spacer)
		mature.setdefault("_".join(genes), set()).add(spacer)

#remove sRNAs that have all spacers in common with another:
mature = {key:val for key, val in mature.items() if val != []}


junctions = ["CTAA", "GACA", "GCAC", "AATC", "GTAA", "TGAA", "ATTA", "CCAG",
 "AGGA", "ACAA", "TAGA", "CGGA", "CATA", "CAGC", "AAGT", "CTCC", "AGAT", "ACCA",
 "AGTG", "GGTA", "GCGA", "AAAA", "ATGA"]
start_junction = "CCTC"
end_junction = "AACG"
repeat_dna = "gtctaagaactttaaataatttctactgttgtagat"
rev_repeat_dna = revcomp(repeat_dna)

if arrays_flag == True:
	oligos = [[] for i in range(int(len(arrays)/len(junctions)+1))]
	out_array_constructs = open("constructs_array.txt",mode="w")
	print("Arrays:")

	previous_gene = ""
	n = 0
	j = 1
	for group in oligos:
		i = 0
		out = open("array_oligo_pool_{0}.txt".format(j),mode="w")
		while i < len(junctions) and n < len(arrays):
			array = arrays[n][1]
			gene = arrays[n][0]
			if gene != previous_gene:
				print("\n")
				print(gene)
			junction = junctions[i]
			revc_junction = revcomp(junction)
			f1 = start_junction + repeat_dna + array[0] + junction
			r1 = revcomp(array[0]) + rev_repeat_dna
			f2 = repeat_dna + array[1]
			r2 = revc_junction + revcomp(array[1]) + rev_repeat_dna + revc_junction
			f3 = junction + repeat_dna + array[2]
			r3 = revcomp(end_junction) + revcomp(array[2]) + rev_repeat_dna
			group.extend([f1,r1,f2,r2,f3,r3])
			print("\n".join([f1,r1,f2,r2,f3,r3]))
			i += 1
			n += 1
			out.write("\n".join([f1,r1,f2,r2,f3,r3]))
			out.write("\n")
			construct="".join([f1,f2,f3,end_junction])
			line = ">" + gene + "_array" + "\n"
			out_array_constructs.write(line)
			line = construct + "\n"
			out_array_constructs.write(line)
		j += 1
		print("--------------------------------------------------------------------------------------\n")

print("Mature gRNAs:\n")
out_mature_constructs = open("constructs_mature.txt",mode="w")

out = open("mature_oligo_pool2.txt",mode="w")
for gene, spacers in mature.items():
	print(gene)
	i = 1
	for spacer in spacers:
		oligo = left_overh + spacer + right_overh
		print(oligo)
		out.write(oligo)
		out.write("\n")
		line = ">" + gene + "_mature_" + str(i) + "\n"
		out_mature_constructs.write(line)
		line = oligo + "\n"
		out_mature_constructs.write(line)
		i += 1
	print("\n")

if non_targeting_flag:	
	print("\nMature gRNAs (non-targeting):")
	out = open("nontargeting_oligo_pool2.txt",mode="w")
	i = 1
	for nt_gRNA in nt_spacers[:min_output_size]:
		oligo = left_overh + nt_gRNA[:spacer_length] + right_overh
		print(oligo)
		out.write(oligo)
		out.write("\n")
		line = "non_targeting_" + str(i) + "\n"
		out_mature_constructs.write(line)
		line = oligo + "\n"
		out_mature_constructs.write(line)
		i += 1

if arrays_flag == True:
	if not any(oligos): #if the list is empty, therefore no arrays have been designed
		print("\nNo arrays could be designed for the input targets.")
	else:
		print("\nNumber of array-cloning reactions:", len(oligos))
		print("\nTarget genes with:\n0 arrays: {0}\n1 array: {1}".format(noarray, onearray))
print("\ntot mature: {0}\n".format(total_mature))

print("Done!")
print(datetime.now() - startTime)
