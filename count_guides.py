import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO
from multiprocessing import Pool, Lock
from multiprocessing.managers import SharedMemoryManager
import argparse

parser = argparse.ArgumentParser(description='Count occurrences of gRNA/arrays in a read file.')

parser.add_argument('infile', type=str, help='Input fastq file.')
parser.add_argument('outfile', type=str, help='Output counts file.')
parser.add_argument('--reference', '-r', type=str, help='Input fastq file. (Default: references.fasta)',
					default = "references.fasta" )
parser.add_argument('--processors', '-p', type=int, help='No. of processors. (Default: 1)', default = 1)

args = parser.parse_args()
infile = args.infile
reference = args.reference
outfile = args.outfile
processors = args.processors

def count_refs(fastq_gen):
# counts each entry in references in the fastq file
	global counter_list
	[_, seq, __] = fastq_gen
	i = -1
	for [_,ref,ref_revc] in references:
		i+=1
		if ref in seq:
			lock.acquire()
			counter_list[i] += 1
			lock.release()
			break
		else:
			if ref_revc in seq:
				lock.acquire()
				counter_list[i] += 1
				lock.release()
				break

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

#initialize list of references
references = []
for record in SeqIO.parse(reference,"fasta"):
	references.append([record.id,str(record.seq.upper()),revcomp(str(record.seq.upper()))])

# count gRNAs:
lock = Lock()
with SharedMemoryManager() as smm:
	counter_list = smm.ShareableList([0 for x in references])
	in_handle = gzip.open(infile,"rt")
	pool=Pool(processors)
	fastq_gen = FastqGeneralIterator(in_handle)
	pool.map(count_refs, fastq_gen)
	pool.close()
	pool.join()
	in_handle.close()
	smm.shutdown()

#write results to outfile:
with open(outfile, mode="w") as outf:
	for ref,count in zip(references,list(counter_list)):
		print(ref[0],count)
		outf.write("\t".join([ref[0],str(count)+"\n"]))