from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import argparse
import sys

parser = argparse.ArgumentParser(description='Extracts gene sequences from gff and genome files and adds the promoter sequence.')

parser.add_argument('gff_file', type=str,
					help='Input annotation (gff) file with sequences of the target genes.')
					
parser.add_argument('genome_file', type=str,
					help='Input genome file.')

parser.add_argument('--outname', '-o', type=str, default="target_genes",
					help='Base name of the output files.'
					'(Default: target_genes.)')

parser.add_argument('--promoter', '-p', default=50, type=int,
					help='Length before the gene start to be considered as the promoter. Set to 0 if you want to avoid targeting '
					'the promoter. (Default: 50)')

parser.add_argument('--genetype', '-g', type=str,
					help='Consider only genes of this type (Third column of the gff file).')			

parser.add_argument('--attribute', '-a', type=str,
					help='Consider only genes that have this string in the attributes column (ninth column of the gff file).'
					'E.g. if you want to filter for the attribute "sRNA_type=Intergenic", write '
					'sRNA_type=Intergenic.')

parser.add_argument('--name', '-n', type=str, default="locus_tag",
					help='Take gene name from this attribute of the gff attribute column (ninth column of the gff file).')	
					
args = parser.parse_args()
gff_file = args.gff_file
genome_file = args.genome_file
outname = args.outname
promoter_length = args.promoter
genetype = args.genetype
attribute = args.attribute
gene_name = args.name

def revcomp(sequence):
#reverse complements sequence
	basecomplement = {'A':'T', 'C':'G', 'T':'A', 'G':'C'} 
	letters = list(sequence) 
	letters.reverse() 
	dna = ''
	for base in letters:
		dna += basecomplement[base] 
	return dna

##### extend sRNA coordinates to add promoter_length before the start:
input_file = open(gff_file, mode="r")
genes = {}
for line in input_file:
	linelist = line.split("\t")
	if line.startswith("#") or not line.rstrip():
		continue
	if genetype is not None:
		if linelist[2] != genetype:
			continue
	if attribute is not None:
		if attribute not in linelist[8]:
			continue
	seqname = linelist[0]
	attributes = {}
	for attr in linelist[8].split(";"):
		attr = attr.split("=")
		attributes[attr[0]] = attr[1]
	try:
		gene = attributes[gene_name].rstrip()
	except KeyError:
		print("ERROR: this line of the gff file: \n\n{}\ndoesn't have the "
		"desired attribute '{}' which is set by the -n flag.".format(line,gene_name))
		sys.exit()
	strand = linelist[6]
	if strand == "+":
		start = str(int(linelist[3]) - promoter_length)
		end = str(int(linelist[4]))
	elif strand == "-":
		end = str(int(linelist[4]) + promoter_length)
		start = str(int(linelist[3]))
	genes[gene] = [seqname,start,end,strand]

with open("{}_coordinates.txt".format(outname), mode="w") as outf:
	for gene,[scaffold,start,end,strand] in genes.items():
		line = "\t".join([scaffold,gene,start,end,strand])
		outf.write(line)
		outf.write("\n")

##### make a fasta file with gene name and sequence:
input_file = open(genome_file)
#make a dictionary out of the genome:
assembly_seqs = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))
#make a list of SeqRecord items of each gene:
genes_sr = []
for gene,coord in genes.items():
	[seqname,start,end,strand] = coord
	gene_sequence = assembly_seqs[seqname][int(start)-1:int(end)].seq.upper()
	if strand == "-":
		gene_sequence = revcomp(gene_sequence)
	gene_sr = SeqRecord(gene_sequence, id=gene, description ="")
	genes_sr.append(gene_sr)
	genes[gene].append(gene_sequence)

#make a dictionary of sequence:geneID for each gene.
#this is to capture genes that have the same sequence
#and store their IDs as a list in the dictionary value
seq_dict = defaultdict(list)
for gene in genes_sr:
	seq_dict[str(gene.seq)].append(gene.id)

#write the output fasta file:
out_fasta = open("{}_sequences.fasta".format(outname), 'w')
for seq, IDs in seq_dict.items():
	#Join the IDs of genes that have the same sequence as gene1|gene2:
	ID = '|'.join(IDs)
	out_fasta.write(">{}\n".format(ID))
	out_fasta.write(seq + "\n")