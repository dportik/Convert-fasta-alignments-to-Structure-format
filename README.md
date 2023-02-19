# Convert fasta alignments to Structure format

This is a script for a specific use-case in which phased alignments have been produced for many genes/loci. The goal is to produce a Structure input file from these alignments by sampling one SNP site per alignment file. These data could come from sequence-capture experiments (UCEs, exons, etc.), but this should also work with ddRADseq loci.

This program requires a set of phased alignments contained in multiple fasta files. Each fasta file should represent a different gene/locus, and each fasta file should have two sequences per individual, representing the phased alleles. Because these data should be phased, there should not be any ambiguity codes present in the sequences. For this script, it does not matter how the individuals are labeled (e.g., `Sample1_a`, `Sample1_b`; `Sample1_0`, `Sample1_1`), though this may matter for Structure downstream. 

This script will identify all the variable columns present in each alignment fasta, and select one column randomly. The final set of sampled alignment columns is equivalent to selecting one SNP per gene/locus. This should ensure linkage is not an issue for the model-based clustering approach in Structure.

The subsampled alignment columns are then translated into the format used by the program Structure. Two outputs are written:

+ `Structure_File_Nucleotides.txt`: contains nucleotide characters `A`, `C`, `T`, `G` (with missing data represented by a `-`).

+ `Structure_File_IntegerCodes.txt`: contains integers `0`, `1`, `2`, `3` (with missing data represented by a `-9`). 

The script is written in Python and is compatible with Python 2 & Python 3.

To use the script:
```
python Convert_Fasta_Set_to_Structure.py -i <path to directory with fasta files> -o <path to output directory>
```

The optional `--remove_singletons` flag can be included to eliminate alignment columns with singleton SNPs from being selected. Those columns are considered "variable" sites but may not contain useful population level information. 
 
Example data are provided in the `Example-Data/Inputs` folder, and will produce the output files contained in the `Example-Data/Outputs` folder. 