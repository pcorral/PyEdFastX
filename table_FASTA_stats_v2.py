#
# Produce a table with statistics from a FASTA & QUAL file
# still incomplete...
#
# Command to execute:
#   $~> python table_FASTA_stats_v#.py your_file.fasta
#
# Essentially, all the preprocessint deals with the fact that what we are
# trying to do is align all sequences by its middle base, so that,
# both extrems in each sequence can be compared
#
#
# Note:
#   1)place this script and get_FASTA_stats.py in the same working.direct
#
#
# Pau C. (initiated 01/04/2012)
##################################################


import sys                     # collect arguments in command line

#from scipy.stats import scoreatpercentile      # use quartiles
#from Bio.SeqIO.QualityIO import PairedFastaQualIterator  # FASTA & QUAL parsers
from numpy import zeros, add            # use statistic functions as
                               # sum, min, max, mean and
                               #array creation and manipulation


#When we import get_FASTA_stats.py, a function [basic_stats()] is executed, and
#two (THREE, update!) values are loaded into the workspace:
#        ->length of the longest sequence
#        ->number of sequences

from get_FASTA_stats_v2 import *
###################################################################

# First we execute the function basic_stats() from the get_FASTA_stats.py file
# to get the two (THREE, update!) variables holding
#        longest seq in the FASTA     -> longest_out
#        number of seq's in the FASTA -> total_seq_out
#        list of sequence length      -> seq_len_list
(longest_out, total_seq_out, seq_len_list) = basic_stats(fasta_file)

# Initialize vectors to hold bases count
A_count_partial = zeros(longest_out, int)
C_count_partial = zeros(longest_out, int)
G_count_partial = zeros(longest_out, int)
T_count_partial = zeros(longest_out, int)
N_count_partial = zeros(longest_out, int)

# Open the FASTA file. Rise an error if not possible.
try:
    handle = open(fasta_file, 'rU')
except IOError:
    print 'cannot open', fasta_file,
    '\nMake sure the file is in the working directory and you \
    have typed it correctly'

# Iterate over only the FASTA file with the Biopython parser in SeqIO.
# We center all sequences in its middle base, being
# len(my_seq)/2 the middle position
# note that odd values divided by 2 return just the integer part as in:
# 5/2 = 2 (2 comes from 2.5)

middle_longest = longest_out / 2  # middle position of longest seq 

for record in SeqIO.parse(handle, "fasta"):

        move_position = -1    # Increment +1 in each iteration
                              # set to -1 to start the iteration with value 0
                              # (0 will point to first position)

        # len(record.seq) evaluates length of each sequence
        position = middle_longest - len(record.seq) / 2

        #Here we iterate over each base in the sequence and add 1
        for base in record.seq:
                move_position += 1
                if base == 'A':
                        A_count_partial[position + move_position] += 1
                elif base == 'C':
                        C_count_partial[position + move_position] += 1
                elif base == 'G':
                        G_count_partial[position + move_position] += 1
                elif base == 'T':
                        T_count_partial[position + move_position] += 1
                elif base == 'N':
                        N_count_partial[position + move_position] += 1
                #here we capture bases that are not an A, C, G, T or N
                #else:
                #    print "base", base, "found\n"
                
                
# add() function from numpy is used to calculate
# total number of bases in each position
calc_matrix = add(A_count_partial, C_count_partial)
calc_matrix = add(calc_matrix, G_count_partial)
calc_matrix = add(calc_matrix, T_count_partial)
per_position_total_bases = add(calc_matrix, N_count_partial)

handle.close()

print per_position_total_bases

# What I have so far is: (which has to be writen in the table [tabulated file])
# Nomenclature -->
#   in_table_name : in_my_code_name = explanation of what it is/ what it holds
#
#    column   : extremes                 = range of base positions
#                                          rang(0:longest_out)
#    count    : per_position_total_bases = bases found in each position
#    A_Count  : A_count_partial          = A base counts in each position
#    C_Count  : C_count_partial          = C base counts in each position
#    G_Count  : G_count_partial          = G base counts in each position
#    T_Count  : T_count_partial          = T base counts in each position
#    N_Count  : N_count_partial          = N base counts in each position
#    Max_count : total_seq_out           = number of sequences
#
# still pending ....

# They only work with the QUAL values:
#   min : min_col                    = Lowest quality score value found in this column.
#   max : max_col                    = Highest quality score value found in this column.
#   sum : sum_col                    = Sum of quality score values for this column.
#   mean : mean_col                  = Mean quality score value for this column.
#   Q1  : Q1_col                     = 1st quartile quality score.
#   med : med_col                    = Median quality score.
#   Q3  : Q3_col                     = 3rd quartile quality score.
#   IQR : IQR_col                    = Inter-Quartile range (Q3-Q1).
#   lW  : lw_col                     = -2*sd. 'Left-Whisker' value (for boxplotting).
#   rW  : rw_col                     = 2*sd. 'Right-Whisker' value (for boxplotting).

#Open the QUAL file
#Note that 'your_file.fasta' is used to get the name of the QUAL file
#Therefore, the FASTA and QUAL file have to have the same name

qual_file_tmp = fasta_file.split(".")[0]         #get the name of the fasta file
qual_file = '.'.join([qual_file_tmp,"qual"])     #join two words by the separator
                                                 #specified 

# Open the QUAL file. Rise an error if not possible.
try:
    handle = open(qual_file, 'rU')
except IOError:
    print 'cannot open', qual_file,
    '\nMake sure the FASTA and QUAL files have the same name and\
    you have typed it correctly'
handle.close()

# Iterate over the QUAL file
###############################

min_col = zeros(longest_out, int)
max_col = zeros(longest_out, int)  
sum_col = zeros(longest_out, int)  
mean_col = zeros(longest_out, int)
Q1_col = zeros(longest_out, int) 
med_col = zeros(longest_out, int)
Q3_col = zeros(longest_out, int) 
IQR_col = zeros(longest_out, int)
lw_col = zeros(longest_out, int) 
rw_col = zeros(longest_out, int)  

tmpL_col = []   # Holds, in order and reset to an empty list when processed,
                # column [base position] of qual values where, 
                # seqs are aligned in its middle position

total_colmn = range(longest_out, 0, -1)  #this produces a list like 
                                        #[10,9,8,7,6,5,4,3,2,1]
                                        #if longest_out=10
iter_num_asc = -1    # keep track of the column 

# Iterate over the positions/columns of the longest seq
for colmn in total_colmn:

    iter_num_asc += 1
    
    # list comprehension is used to check the list of seq's length
    # and see which positions fulfill the length given by colmn
    # if fulfills, the seq has a base in this column
    
    index_fulfill = [i for i,x in enumerate(seq_len_list) if x >= colmn]
    index_fulfill.sort()    # in ascending order. Inplace replacement
    
    print 'those seqs that fulfill', index_fulfill
    # 1)keep passing entries in the QUAL file, 
    # until an index matches with the entry
    
    # 2)The file is opened, parsed from top to bottom
    # and closed just after the sequence is found
    for pos_val in index_fulfill:
        handle = open(qual_file, 'rU')
        next_entry = -1  #it will increment in every iteration
                     #when next_entry=0 it will point to the first entry
    
        # Here we reach the sequence of QUAL values and calculate some stats
        # that are saved in the temporary list tmpL_col
        for record in SeqIO.parse(handle, "qual"):
            next_entry += 1
            
            # only manipulate the seqs given by pos_val
            flag_stop_parsing = 0
            if pos_val == next_entry:
                flag_stop_parsing = 1    
                # This is the way I found to access the phred qual values 
                # note that only the value of interest is retrieved by indeximg
                print record.letter_annotations['phred_quality'][seq_len_list[pos_val] - colmn]
                print 'value found in this position', seq_len_list[pos_val] - colmn
            if flag_stop_parsing == 1:
                print 'hola4'
                break
                #tmpL_col.append(record.letter_annotations['phred_quality']
                #                                [pos_val - colmn])
                
                
    #print 'min value is', min(tmpL_col)
    handle.close()



    
        
        


    
    
    
    