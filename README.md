# CRISPRi_designer
Tools to design gRNA/array pools targeting a list of genes. 

## Dependencies
Create a conda environment from the environment.yml file:
```
conda env create -f environment.yml
```

## get_sequences_from_gff.py
The first step in the pipeline is to obtain a fasta file with the sequences of the genes that are to be targeted. These sequences should contain the promoter. The file can be created from a gff annotation with the get_sequences_from_gff.py script.

### Usage
```
usage: get_sequences_from_gff.py [-h] [--outname OUTNAME]
                                 [--promoter PROMOTER] [--genetype GENETYPE]
                                 [--attribute ATTRIBUTE] [--name NAME]
                                 gff_file genome_file

Extracts gene sequences from gff and genome files and adds the promoter
sequence.

positional arguments:
  gff_file              Input annotation (gff) file with sequences of the
                        target genes.
  genome_file           Input genome file.

optional arguments:
  -h, --help            show this help message and exit
  --outname OUTNAME, -o OUTNAME
                        Base name of the output files. The .fasta extension
                        will be appended to this. (Default: target_genes.)
  --promoter PROMOTER, -p PROMOTER
                        Length before the gene start to be considered as the
                        promoter. (Default: 50)
  --genetype GENETYPE, -t GENETYPE
                        Consider only genes of this type (Third column of the
                        gff file).
  --attribute ATTRIBUTE, -a ATTRIBUTE
                        Consider only genes that have this attribute (ninth
                        column of the gff file).E.g. if you want to filter for
                        the attribute "sRNA_type=Intergenic", write
                        sRNA_type=Intergenic.
  --name NAME, -n NAME  Take gene name from this attribute of the gff
                        attribute column (ninth column of the gff file).
```
### Input
* A gff annotation (*gff_file*) in the standard format (tab-separated, no header, columns: seqname,feature,start,end,score,strand,frame,attribute).
* A (multi)fasta file (*genome_file*) containing the genome file (chromosome(s) and plasmid(s), if present).

### Output
* A *outname*_coordinates.txt file containing the coordinates of each gene (promoter included).
* A *outname*_sequences.fasta file containing the genes' sequences (promoter included).

*outname* is defined by the --outname option.

## PAM_frequency.py
(Optional). This script scans the input sequences to find the frequency of selected PAM sequences among the input genes.

### Usage

```
usage: PAM_frequency.py [-h] [--input_file INPUT_FILE] [--PAMs PAMS]
                        [--spacer_length SPACER_LENGTH] [--binsize BINSIZE]

Calculates the frequency of each PAM in the input fasta file.

optional arguments:
  -h, --help            show this help message and exit
  --input_file INPUT_FILE, -i INPUT_FILE
                        Fasta file with sequences of the target genes +
                        promoters. (Default: target_genes_sequences.fasta).
  --PAMs PAMS, -pam PAMS
                        List of PAM sequences and their position relative to
                        the protospacer. Format: PAM1-position,PAM2-position
                        (can be more than two). (Default: TTV-5,NGG-3).
  --spacer_length SPACER_LENGTH, -sl SPACER_LENGTH
                        Length of each spacer. No more than 26 if the -a flag
                        is used. (Default: 20).
  --binsize BINSIZE, -b BINSIZE
                        Size of the bins in the final heatmap. (Default: 5).
```
### Input
* A fasta file (*input_file*) file containing the sequences of the target genes (promoter included).

### Output
* *PAM_frequency_counts.tsv*: tab-separated file with one PAM per row and number of genes that have 0...n occurrences of the PAM. Column names are the number of PAMs/gene.
* *PAM_frequency_heatmap.pdf*: Heatmap with the same results, binned in bins of size set by *--binsize*.

## design_CRISPRi_gRNAs.py
This script designs the gRNA/array library. PAM sequence (*--PAM*), strand preference inside the transcribed region (*PAM_preference*) and orientation respective to the protospacer (*--PAM_orientation*) can be freely set. gRNAs in the promoter can target both strands.

For Cas12a, the script can design 1 array per gene if the *--arrays* flag is set. In this case, the oligonucleotides designed for the arrays are structured for cloning via CRATES (PMID: 31270316). The default *--left_overlap* and *--right_overlap* overlap sequences are for cloning into AWP-029.

This script support multiprocessing (*--processors*), however >1 processors are only used when looking for offtargets and designing arrays (*--arrays*). The second case is the one where multiple processors really make a difference. With 40 processors, the script takes about 30 minutes for the ~150 genes of ~150 nt length.


### Usage

```
usage: design_CRISPRi_gRNAs.py [-h] [--input_file INPUT_FILE] [--coordinates COORDINATES] [--reference REFERENCE] [--spacer_length SPACER_LENGTH]
                               [--promoter_length PROMOTER_LENGTH] [--non_targeting [NON_TARGETING]] [--PAM PAM] [--PAM_preference {template,nontemplate}]
                               [--PAM_orientation {5prime,3prime}] [--spacers SPACERS] [--arrays] [--folding FOLDING] [--left_overlap LEFT_OVERLAP]
                               [--right_overlap RIGHT_OVERLAP] [--processors PROCESSORS]

Designs spacers targeting the input genes.

optional arguments:
  -h, --help            show this help message and exit
  --input_file INPUT_FILE, -i INPUT_FILE
                        Fasta file with sequences of the target genes + promoters. (Default: target_genes_sequences.fasta)
  --coordinates COORDINATES, -c COORDINATES
                        File with coordinates of the target genes. (Default: target_genes_coordinates.txt)
  --reference REFERENCE, -r REFERENCE
                        Genome sequence file (Default: genome.fasta).
  --spacer_length SPACER_LENGTH, -sl SPACER_LENGTH
                        Length of each spacer. No more than 26 if the -a flag is used. (Default: 20)
  --promoter_length PROMOTER_LENGTH, -pl PROMOTER_LENGTH
                        Length of the promoter region. (Default: 50)
  --non_targeting [NON_TARGETING], -nt [NON_TARGETING]
                        The script designs nontargeting gRNAs if this flag is present. If only -nt is set, the number of designed nontargeting gRNAs is the largest number
                        between 20 and total number of gRNAs divided by 2. If -nt n is set, the script designs n nontargeting gRNAs.
  --PAM PAM, -pam PAM   PAM sequence. (Default: TTV)
  --PAM_preference {template,nontemplate}, -pp {template,nontemplate}
                        Which strand within the coding region should be targeted? (Default: template)
  --PAM_orientation {5prime,3prime}, -po {5prime,3prime}
                        At which end of the protospacer is the PAM located? (Default: 5prime)
  --spacers SPACERS, -s SPACERS
                        Number of spacers (including one array if using the -a flag) to be designed per targeted gene. (Default: 4)
  --arrays, -a          The script attempts to design one array per target gene if this flag is present. This flag can be used only if --PAM is TTV.
  --folding FOLDING, -f FOLDING
                        Minimal accepted correct folding probability of the arrays. Scale: 0-1. (Default: 0.2)
  --left_overlap LEFT_OVERLAP, -lo LEFT_OVERLAP
                        Left overhang for cloning grna spacers. (Default: atctttgcagtaatttctactgttgtagat)
  --right_overlap RIGHT_OVERLAP, -ro RIGHT_OVERLAP
                        Right overhang for cloning grna spacers. (Default: ccggcttatcggtcagtttcacctgattta)
  --processors PROCESSORS, -p PROCESSORS
                        Number of processors. A real impact of using >1 processors is seen only when also using the -a option. (Default: 1)
```

### Input
The two files generated by _get_sequences_from_gff_ with sequences and coordinates of the genes to be targeted (plus the promoter). if different from the default, the file names can be indicated by *--input_file* and *--coordinates*.

### Output
The following files are generated by the script:
* *references.fasta*: Fasta file with all gRNAs (including nontargeting if the *-nt* option is used) and arrays (if the *-a* option is used).
* *gRNA_oligo_pool.txt*: text file with all gRNAs + homology arms (sequences only). Comprises the non targeting gRNAs too if the *-nt* option is used. Might be useful for ordering the oligonucleotides.
* *array_oligo_pool_n.txt*: only if the *-a* option is used. Text file(s) with the sequences of the oligos needed to construct the array via CRATES. One file is generated per each pool to be ordered.

## count_guides.py
This script makes a count table file with the counts of each gRNA/array in a fastq reads file. The input read file must be one (i.e. a unique file with merged mate reads if sequencing was paired-end).

### Usage

```
usage: count_guides.py [-h] [--reference REFERENCE] [--processors PROCESSORS] infile outfile

Count occurrences of gRNA/arrays in a read file.

positional arguments:
  infile                Input fastq file.
  outfile               Output counts file.

optional arguments:
  -h, --help            show this help message and exit
  --reference REFERENCE, -r REFERENCE
                        Input fastq file. (Default: references.fasta)
  --processors PROCESSORS, -p PROCESSORS
                        No. of processors. (Default: 1)

```
### Input
* A gzipped fastq file (*infile*) containing reads from a sequencing experiment.

### Output
* *outfile*: a tab-separated file with one colum with the name of the gRNA/array (from *reference*) and a second column with its count number.

## Examples

* Extract sequences of the targeted genes and promoters (of length 35 nt). The genome is genome.fasta, the gff file is annotation.gff. We want to filter the gff file for CDS entries only, while the name should be taken from the "gene_name" attribute:

	```
	get_sequences_from_gff.py -p 35 -t CDS -a gene_name annotation.gff genome.fasta 
	```

* Calculate frequencies of the following PAMs: TTV (5 prime of the protospacer), NGG (3 prime), NNNNGATT (3 prime) and NNGRR (3 prime). The target sequences are in targets.fasta and we want the final heatmap results to be binned in bins of size 3:

	```
	PAM_frequency.py -i targets.fasta -pam TTV-5,NGG-3,NNNNGATT-3,NNGRR-3 -b 3
	```

* Design a library of gRNAs, arrays and nontargeting gRNAs for Cas12a, with PAM TTV, template strand preference and 5prime orientation. The genome is genome.fasta and the target genes are found in target_genes_sequences.fasta and target_genes_sequences.fasta:

	```
	design_CRISPRi_gRNAs.py -i target_genes_sequences.fasta -c target_genes_sequences.fasta -r genome.fasta -a -nt -pam TTV -pp template -po 5prime
	```
 
* Design a library of gRNAs for Cas9, with PAM NGG, nontemplate strand preference and 3prime orientation. The genome is genome.fasta and the target genes are found in target_genes_sequences.fasta and target_genes_sequences.fasta:
	
	```
	design_CRISPRi_gRNAs.py -i target_genes_sequences.fasta -c target_genes_sequences.fasta -r genome.fasta -pam NGG -pp nontemplate -po 3prime
	```

* Design a library of targeting and nontargeting gRNAs for Cas9, with PAM NGG, nontemplate strand preference and 3prime orientation. We want 50 nontargeting gRNAs.
	
	```
	design_CRISPRi_gRNAs.py -nt 50 -pam NGG -pp nontemplate -po 3prime
	```

* Count the gRNAs and arrays of the file *references.fasta* in the fastq *reads.fastq* and save to *counts.tsv*. Use 20 processors.

	```
	python count_guides.py -p 20 -r references.fasta reads.fastq counts.tsv
	```