# Convert fasta alignments to Structure format

This is a script for a specific use-case in which phased alignments have been produced for many genes/loci. The goal is to produce a Structure input file from these alignments by sampling one SNP site per alignment.

This program requires a set of phased alignments contained in multiple fasta files. Each fasta file should represent a different gene/locus, and each fasta file should have two sequences per individual, representing the phased alleles. 

This program will identify all the variable columns present in each alignment fasta, and select one column randomly. The final set of sampled alignment columns is equivalent to selecting one SNP per gene/locus. 

The selected alignment columns are then translated into the format used by the program Structure. Two outputs are written, one which contains nucleotide characters (with missing data represented by a `-`), and one which contains integers instead (with missing data represented by a `-9`). 

The script is written in Python and is compatible with Python 2 & Python 3.

To use the script:
```
python Convert_Fasta_Set_to_Structure.py -i <path to directory with fasta files> -o <path to output directory>
```
 
Example data are provided in the `Example-Data/Inputs` folder, and will produce the output files contained in the `Example-Data/Outputs` folder. 