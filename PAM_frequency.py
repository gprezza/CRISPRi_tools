from Bio import SeqIO
from Bio import Data
from itertools import product
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('agg')
from collections import OrderedDict
import argparse

parser = argparse.ArgumentParser(description='Calculates the frequency of each PAM in the input fasta file.')

parser.add_argument('--input_file', '-i', type=str, default="target_genes_sequences.fasta",
					help='Fasta file with sequences of the target genes + promoters. (Default: target_genes_sequences.fasta).')

parser.add_argument('--PAMs', '-pam', default="TTV-5,NGG-3", type=str, help='List of PAM sequences and their position '
					'relative to the protospacer. Format: PAM1-position,PAM2-position (can be more than two). (Default: TTV-5,NGG-3).')

parser.add_argument('--spacer_length', '-sl', default="20", type=int,
					help='Length of each spacer. No more than 26 if the -a flag is used. (Default: 20).')

parser.add_argument('--binsize', '-b', default="5", type=int,
					help='Size of the bins in the final heatmap. (Default: 5).')

args = parser.parse_args()
input_file = args.input_file
PAM_args = args.PAMs
spacer_length = args.spacer_length
binsize = args.binsize

ambiguous_alph = {"N":["A","C","G","T"], "V":["A","C","G"], "H":["A","C","T"], "D":["A","G","T"],
	"B":["C","G","T"], "M":["A","C"], "K":["G","T"], "W":["A","T"], "S":["C","G"], "Y":["C","T"],
	"R":["A","G"], "A": ["A"], "C": ["C"], "G": ["G"], "T": ["T"]}

def remove_ambiguous(sequence):
# returns a list of all possible sequences
# given the ambiguous sequence input (e.g. with N,V, etc. bases)
	return ["".join(i) for i in product(*[ ambiguous_alph[j] for j in sequence ]) ]

def revcomp(sequence):
#reverse complements sequence
	complement = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N','Y':'R',
	'R':'Y', 'W':'W', 'S':'S', 'K':'M', 'M':'K', 'D':'H', 'H':'D',
	'V':'B', 'B':'V'} 
	bases = list(sequence) 
	bases.reverse() 
	revcomp = ''
	for base in bases:
		revcomp += complement[base] 
	return revcomp

def find_all(string, substr):
# makes a list with the indexes of all occurrences of substr in string
	listindex=[]
	start = 0
	i = string.find(substr, start)
	while i >= 0:
		listindex.append(i)
		i = string.find(substr, i + 1)
	return listindex	

def find_guides(seq):
# finds gRNAs in seq
	guides_bottom = []
	guides_top = []
	for PAM_sequence in PAM_list:
		guides_top.extend(find_all(seq, revcomp(PAM_sequence)))
		guides_bottom.extend(find_all(seq, PAM_sequence))
	return guides_top, guides_bottom

def binned_sum(data):
# Sums the values in data in bins of binsize size. Returns an OrderedDict
# where keys are the bin number (starting from 1) and values the sum of the
# values in that bin. The first item (key=0) equals to data[0].
	binned_data = OrderedDict()
	binned_data[0] = data[0]
	j = 0
	for i in bins:
		binned_data[j+1] = 0
		for k in range(bins[j]-binsize+1,bins[j]+1):
			if k in data:
				binned_data[j+1] += data[k]
		j += 1
		if j == len(bins):
			return binned_data
	return binned_data

#create a list of the PAMs in the style [["TTV",5],["NGG",3]]

PAM_args = PAM_args.split(",")
PAMs = []
for PAM in PAM_args:
	seq,pos = PAM.split("-")
	pos = int(pos)
	if pos not in (3, 5):
		print("ERROR: The PAM position you gave to the --PAMs flag is invalid. It has to be either 3 or 5.")
		exit()
	PAMs.append([seq,pos])

#make a PAMs_dic dictionary with each PAM as key and as value a dictionary with 
#key: number of gRNAs/gene value:number of gRNAs that have that number
PAMs_dic = {}
max_no = 0
for PAM in PAMs:
	spacer_count_dict = OrderedDict()
	sRNA_list = SeqIO.parse(input_file, 'fasta')
	PAM_list = remove_ambiguous(PAM[0])
	PAM_length = len(PAM_list[0])
	for sRNA in sRNA_list:
		sRNA_sequence = str(sRNA.seq)
		sRNA_length = len(sRNA_sequence)
		guides_top, guides_bottom = find_guides(sRNA_sequence)
		gRNA_number = 0
		if PAM[1] == 3: #eg NGG
			for pos in guides_top:
				if pos <= sRNA_length - spacer_length - PAM_length:
					gRNA_number += 1
			for pos in guides_bottom:
				if pos >= spacer_length:
					gRNA_number += 1
		else: #eg TTV
			for pos in guides_bottom:
				if pos <= sRNA_length - spacer_length - PAM_length:
					gRNA_number += 1
			for pos in guides_top:
				if pos >= spacer_length:
					gRNA_number += 1
		if gRNA_number in spacer_count_dict:
			spacer_count_dict[gRNA_number] += 1
		else:
			spacer_count_dict[gRNA_number] = 1
	PAMs_dic[PAM[0]] = spacer_count_dict
	max_no = max(max_no,max(spacer_count_dict))

#add an entry of 0 for each value that is missing
for PAM, counts in PAMs_dic.copy().items():
	for i in range(0,max_no+1):
		if i not in counts:
			PAMs_dic[PAM][i] = 0

#Create the the results table file
outfile = open("PAM_frequency_counts.tsv",mode="w+")
outfile.write("\t".join(["PAM"] + [str(x) for x in range(max_no+1)]))
outfile.write("\n")

#Bin the dictionaries in PAMs_dic into bins of binsize
#at the same time, write to file
num_bins = max_no // binsize + 1
bins = [i * binsize for i in range(1,num_bins+1)]
max_no_binned = 0
for PAM, counts in PAMs_dic.copy().items():
	ordered_counts = OrderedDict(sorted(counts.items(),reverse=False))
	outfile.write("\t".join([PAM] + [str(x) for x in ordered_counts.values()]))
	outfile.write("\n")
	PAMs_dic[PAM]=binned_sum(ordered_counts)
	max_no_binned = max(max_no_binned,max(PAMs_dic[PAM].values()))

#prepare lists for the heatmap
x_values = [] #x values in the final heatmap
y_values = [] #y values in the final heatmap
y_labels = [] #y_labels in the final heatmap
x_labels = [str(0)] #x_labels in the final heatmap

for PAM in PAMs:
	pam = PAM[0]
	x_values.append(list(PAMs_dic[pam].keys()))
	y_values.append(PAMs_dic[pam].values())
	y_labels.append(pam)
y_values=[list(y) for y in y_values]

for i in bins:
	x_labels.extend(["{}-{}".format(i-binsize+1,i)])

# Create the heatmap
plt.rcParams.update({"font.size": 13, "xtick.labelsize": 13, "ytick.labelsize": 13, "axes.labelsize": 13})
plt.figure(figsize=(len(x_labels), len(y_labels)))
y = np.array(y_values)
heatmap = plt.pcolormesh(y, cmap=plt.cm.Blues, vmax=max_no_binned)

# Add annotations
for i in range(len(y_labels)):
    for j in range(len(x_labels)):
        plt.text(j + 0.5, i + 0.5, str(y[i, j]), ha='center', va='center')

# Customize, ticks, labels and titles
plt.yticks(np.arange(len(y_labels)) + 0.5, y_labels)
plt.xticks(np.arange(len(x_labels)) + 0.5, x_labels, rotation=30, ha='right',rotation_mode='anchor')
plt.colorbar(label='Number of gRNAs')
plt.xlabel("Number of PAMs/gene")
plt.gca().invert_yaxis()

# Adjust the plot limits
plt.xlim(0, len(x_labels))
plt.ylim(0, len(y_labels))

# Add boundary lines
plt.axhline(y=0, color='0.2', linewidth=2)
plt.axhline(y=len(y_labels), color='0.2', linewidth=2)
plt.axvline(x=0, color='0.2', linewidth=2)
plt.axvline(x=len(x_labels), color='0.2', linewidth=2)

plt.savefig('PAM_frequency_heatmap.pdf', bbox_inches='tight')