# CRISPRi_designer
Tools to design gRNA/array pools targeting a list of genes. 

## Dependencies
Create a conda environment from the environment.yml file:
```
conda env create -f environment.yml
```
## get_sequences_from_gff
The first step in the pipeline is to obtain a fasta file with the sequences of the genes that are to be targeted. These sequences must contain the promoter. The file can be created from a gff annotation with the get_sequences_from_gff.py script.

### Usage
```
usage: get_sequences_from_gff.py [-h] [--outname OUTNAME]
                                 [--promoter PROMOTER] [--genetype GENETYPE]
                                 [--attribute ATTRIBUTE] [--name NAME]
                                 gff_file fasta_file

Extract gene sequences from gff and genome files.

positional arguments:
  gff_file              Input annotation (gff) file with sequences of the
                        target genes.
  fasta_file            Input genome file.

optional arguments:
  -h, --help            show this help message and exit
  --outname OUTNAME, -o OUTNAME
                        Base name of the output files. The .fasta extension
                        will be appended to this. (Default: target_genes.
  --promoter PROMOTER, -p PROMOTER
                        Length before the gene start to be considered as the
                        promoter. (Default: 50)
  --genetype GENETYPE, -t GENETYPE
                        Consider only genes of this type (Third column of the
                        gff file).
  --attribute ATTRIBUTE, -a ATTRIBUTE
                        Consider only genes that have this attribute (ninth
                        column of the gff file).E.g. if you want to filter for
                        the attribute "sRNA_type=Intergenic"
                        writesRNA_type=Intergenic.
  --name NAME, -n NAME  Take gene name from this attribute of the gff
                        attribute column (ninth column of the gff file)
```
### Input files
* A gff annotation (*gff_file*) in the standard format (tab-separated, no header, columns: seqname,feature,start,end,score,strand,frame,attribute).
* A (multi)fasta file (*fasta_file*) containing the genome file (chromosome(s) and plasmid(s), if present).

### Output
* A _basename_coordinates.txt_ file containing the coordinates of each gene (promoter included).
* A _basename_sequences.fasta_ file containing the genes (promoter included).
_basename_ is defined by the --outname option.

## design_CRISPRi_gRNAs.py
This script designs the gRNA/array library.

### Usage

```
usage: design_CRISPRi_gRNAs.py [-h] [--processors PROCESSORS]
                               [--in_name IN_NAME] [--coordinates COORDINATES]
                               [--folding FOLDING] [--left LEFT]
                               [--right RIGHT] [--arrays] [--non_targeting]
                               [--PAM PAM] [--spacer_length SPACER_LENGTH]
                               [--promoter_length PROMOTER_LENGTH]
                               [--mutants MUTANTS] --reference REFERENCE

Designs spacers targeting the input genes.

optional arguments:
  -h, --help            show this help message and exit
  --processors PROCESSORS, -p PROCESSORS
                        Number of processors (Default: 1).
  --in_name IN_NAME, -i IN_NAME
                        Input base name of the fasta file with sequences of
                        the target genes.(Default: target_genes.)
  --coordinates COORDINATES, -c COORDINATES
                        File with coordinates of the target genes. (Default: "
                        --in_name + _coordinates.txt".)
  --folding FOLDING, -f FOLDING
                        Minimal accepted correct folding probability of the
                        arrays. Scale: 0-1. (Default: 0.2).
  --left LEFT, -lo LEFT
                        Left overhang for cloning mature spacers. (Default:
                        atctttgcagtaatttctactgttgtagat).
  --right RIGHT, -ro RIGHT
                        Right overhang for cloning mature spacers. (Default:
                        ccggcttatcggtcagtttcacctgattta).
  --arrays, -a          The script attempts to design one array per target
                        gene if this flag is present.
  --non_targeting, -nt  The script designs nontargeting gRNAs if this flag is
                        present.
  --PAM PAM, -pam PAM   PAM sequence. (Default: TTV).
  --spacer_length SPACER_LENGTH, -sl SPACER_LENGTH
                        Length of each spacer. No more than 26 if the -a flag
                        is used. (Default: 20).
  --promoter_length PROMOTER_LENGTH, -pl PROMOTER_LENGTH
                        Length of the promoter region. (Default: 50).
  --mutants MUTANTS, -m MUTANTS
                        Number of spacers (including one array if using the -a
                        flag) to be designed per targeted gene. (Default: 4).
  --reference REFERENCE, -r REFERENCE
                        Folder containing the genome sequence(s).
```
### Input
* The two files generated by _get_sequences_from_gff_, which are indicated by *in_name* (and *--coordinates*, if the names are different).

### Output
The following files are generated by the script:
* _constructs_mature.txt_: fasta file with all gRNAs + homology arms. Non-targeting gRNAs are also included if the _-nt_ option is used.
* _mature_oligo_pool2.txt_: text file with all gRNAs + homology arms. Like the previous file, but without sequence ID.
* _nontargeting_oligo_pool2.txt_: only if the _-nt_ option is used. Text file with the non-targeting gRNAs + homology arms.
* _constructs_array.txt_: only if the _-a_ option is used. Fasta file with all arrays.
* _array_oligo_pool_n.txt_: only if the _-a_ option is used. n Text files with the sequences of the oligos needed to construct the array via CRATES. One file is generated per each pool to be ordered.