import argparse
import os
import shutil
import random
import operator
from datetime import datetime

def get_args():
    '''
    Get arguments from command line.
    '''
    parser = argparse.ArgumentParser(
            description="""Convert_Fasta_Set_to_Structure - A script to convert""")
    parser.add_argument("-i", "--in_dir",
                            required=True,
                            help="REQUIRED: The full path to a directory which contains the "
                            "input fasta files.")
    parser.add_argument("-o", "--out_dir",
                            required=True,
                            help="REQUIRED: The full path to an existing directory to write "
                            "output files.")
    return parser.parse_args()        

def gather_taxa(f_list):
    """
    Iterate over fasta files in f_list to obtain all 
    unique sample/taxon names in the file set. Returns 
    a sorted list of these names.
    """
    # initiate set to store unique sample/taxon names
    ftaxon_set = set()
    for f in f_list:
        # get sample/taxon names from fasta
        fset = specific_taxon_finder(f)
        # add to set
        ftaxon_set.update(fset)
    # turn set into sorted list
    ftaxon_list = sorted(ftaxon_set)
    print("\n\nFound {} taxon names across {} fasta files.\n\n".format(len(ftaxon_list)/2, len(f_list)))
    
    return ftaxon_list

def specific_taxon_finder(f):
    """
    Function to obtain all the fasta description lines 
    in fasta file f. These descriptions should be a single 
    string with no spaces. Returns the set of names. 
    """
    with open(f, 'r') as fh:
        full_taxa = set([l.replace(">",'').replace('\n','') for l in fh if l.startswith(">")])
        
    return full_taxa

def collect_dicts(f_list, taxon_list, select):
    """
    Identifies variable columns in each fasta file and selects 
    one column randomly to obtain a SNP site for the fasta. Creates 
    a new dictionary for each fasta with sample/taxon names as keys 
    and column base as values. Returns a list of these dictionaries. 
    """
    dict_list, alns_used = [], []
    # iterate over fasta file list
    for f in f_list:
        # convert fasta into dictionary structure
        f_dict = fasta_dict(f, taxon_list)
        # get a list of all sequences sorted by order of sample/taxon
        aln = collect_aln(f_dict, taxon_list)
        # select a variable column from the alignment
        # if one found will be a list of the bases present,
        # else will return None
        snp = process_aln(aln, select)
        # if variable column found
        if snp is not None:
            # create new dictionary with sample/taxon as key and
            # SNP site as value (from the aln column selected)
            snp_dict = new_dict(snp, f_dict)
            # append this dictionary to a list
            dict_list.append(snp_dict)
            # append fasta name to list to track which fastas are used
            alns_used.append(f)
        # if no variable columns found, have to skip fasta
        elif snp is None:
            print("\tSkipping locus!")
            
    return dict_list, alns_used

def fasta_dict(f, taxon_list):
    """
    Function to convert any fasta file (f) into
    dictionary structure with taxon as key and 
    sequence as value. Operates using a taxon_list 
    which indicates the set of expected samples/taxa. 
    If not present, sample/taxon receives missing data 
    value instead of actual sequence. Returns dictionary. 
    """
    print("\nParsing file {}:".format(f))
    aln = []
    f_dict = {}
    with open(f, 'r') as fh:
        lines = [l.strip() for l in fh if l.strip()]
    for line in lines:
        if line.startswith(">"):
            new_key = line.replace(">",'').replace('\n','')
            f_dict[new_key] = ""
        else:
            f_dict[new_key] += line.upper()

    seqlen = len(f_dict[random.choice(list(f_dict.keys()))])
    missing = seqlen * '-'
    
    for t in taxon_list:
        if t in f_dict:
            pass
        else:
            f_dict[t] = missing
            
    return f_dict

def collect_aln(f_dict, taxon_list):
    """
    Obtain a list of the sequences present in this 
    fasta file based on the fasta dictionary. Seqs will 
    be in order sorted by sample/taxon name. Returns 
    list of sequences. 
    """
    aln = []
    for i in sorted(f_dict):
        #print i, f_dict[i]
        aln.append(f_dict[i])
        
    return aln

def process_aln(aln, select):
    """
    Takes list of sequences (aln) produced by collect_aln(). 
    Checks each column to see if they contain 2 or more bases 
    which indicates variable site. If so, stores column in a list. 
    If at least one variable column, selects one randomly to process 
    further, and returns list of bases in column. If no variable 
    columns are found, returns None instead. 
    """
    kept_cols = []
    for i in range(0, len(aln[0])):
        column = []
        for a in aln:
            #print a[i]
            column.append(a[i])
        # determine number of unique & valid nucleotide characters in column 
        flattened, miabp = summarize_col(column)
        if len(flattened) >= 2:
            kept_cols.append([i, miabp])
            #print "Column {} is variable and has {} missing bases".format(i, miabp)
    if kept_cols:
        print("\tFound {} variable sites to choose from.".format(len(kept_cols)))
        if select == "random":
            SNP = random.choice(kept_cols)
        elif select == "sort":
            kept_cols.sort(key=operator.itemgetter(1))
            SNP = kept_cols[0]
        #print SNP
        print("\tSelected alignment column {}.".format(int(SNP[0])+1))
        site = SNP[0]
    else:
        print("\tFound {} variable sites to choose from.".format(len(kept_cols)))
        site = None
        
    return site

def summarize_col(col):
    """
    Function to clean an alignment column (col). Removes 
    missing data characters and counts missing data. Returns 
    a set of the characters in the column (minus missing data 
    characters) and a count of missing data characters. 
    """
    flattened = set(col)
    bad_chars = ["N", "-", "?"]
    for b in bad_chars:
        if b in flattened:
            flattened.remove(b)
    ncount = sum(1 for c in col if c == "N")
    dcount = sum(1 for c in col if c == "-")
    qcount = sum(1 for c in col if c == "?")
    miabp = int(0)
    miabp += ncount
    miabp += dcount
    miabp += qcount

    return flattened, miabp

def new_dict(snp, f_dict):
    """
    Create dictionary structure from a list of the 
    bases in a variable alignment column and the fasta 
    dictionary. New dictionary (snp_dict) has sample/taxon 
    as key and base from the chosen alignment column as value. 
    Returns the new dictionary.
    """
    snp_dict = {}
    for i in f_dict:
        snp_dict[i] = f_dict[i][snp]
        #print f_dict[i][snp]
        
    return snp_dict
        
def concatenate(dict_list, taxon_list):
    """
    This function creates a new dictionary where the 
    keys are sample/taxon names and the values are a 
    concatenated string of all the individual bases selected 
    from the sampled alignment columns from the fasta files. 
    Returns the new dictionary.
    """
    print("\nConverting format of SNP data collected for all taxa...\n")
    concat_dict = {}
    for t in taxon_list:
        concat_dict[t] = ""
        for d in dict_list:
            concat_dict[t] += " {}".format(d[t])
    return concat_dict
        
def convert_format(concat_dict, taxon_list):
    print("\nGathering SNP data collected for all taxa...\n")
    final_dict = {}
    for t in taxon_list:
        newvals = concat_dict[t].replace('A','0').replace('C','1').replace('T','2').replace('G','3').replace('-','-9').replace('N','-9').replace('?','-9')
        final_dict[t] = newvals
    return final_dict

def write_output(d, taxon_list, label):
    with open(label, 'a') as fh:
        for t in taxon_list:
            fh.write("{}{}\n".format(t,d[t]))

def main():
    tb = datetime.now()
    args = get_args()
    os.chdir(args.in_dir)
    # obtain a list of the input fasta files
    f_list = sorted([f for f in os.listdir('.') if f.endswith(('.fa', '.fna', '.fasta'))])
    # get list of samples/taxa found in the fasta files
    taxon_list = gather_taxa(f_list)
    # sample one variable alignment column per fasta
    dict_list, alns_used = collect_dicts(f_list, taxon_list, "random")
    # obtain a concatenated sequence of the sampled bases from the chosen columns
    concat_dict = concatenate(dict_list, taxon_list)
    # convert the nucleotide characters into structure appropriate codes (numerical)
    final_dict = convert_format(concat_dict, taxon_list)
    print("\nSuccessfully obtained a SNP site from {} of {} fasta files.\n".format(len(alns_used), len(f_list)))

    # move to output directory and write
    os.chdir(args.out_dir)
    write_output(concat_dict, taxon_list, "Structure_File_Nucleotides.txt")
    write_output(final_dict, taxon_list, "Structure_File_IntegerCodes.txt")
    print("\nWriting output files.\n")
    tf = datetime.now()
    print("\n\nFinished. Total elapsed time: {0} (H:M:S)\n\n".format(tf - tb))
    
    
if __name__ == '__main__':
    main()
