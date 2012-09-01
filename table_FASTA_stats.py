# Produce a table with statistics from a FASTA & QUAL file in the form
# suitable for FastX_toolkit [http://hannonlab.cshl.edu/fastx_toolkt]
#
# How the information is managed in the scrip:
# After aligning all sequences by its middle base, the following statistics
# are calculated:
# Nomenclature -->
#       Name in      Name in       explanation of what it is/
#   the tab-table :  my code   =        what it holds
#   -----------------------------------------------------------
#
#   column   : range(0, max_seq_len) = range of base positions
#   count    : pos_total_bases = bases found in each position
#   A_Count  : A_count_partial = A base counts in each position
#   C_Count  : C_count_partial = C base counts in each position
#   G_Count  : G_count_partial = G base counts in each position
#   T_Count  : T_count_partial = T base counts in each position
#   N_Count  : N_count_partial = N base counts in each position
#   Max_count : total_seq_out  = number of sequences in total in the file
#   min  : min_col  = Lowest quality score value found in this column.
#   max  : max_col  = Highest quality score value found in this column.
#   sum  : sum_col  = Sum of quality score values for this column.
#   mean : mean_col = Mean quality score value for this column.
#   Q1   : Q1_col   = 1st quartile quality score.
#   med  : med_col  = Median quality score.
#   Q3   : Q3_col   = 3rd quartile quality score.
#   IQR  : IQR_col  = Inter-Quartile range (Q3-Q1).
#   lW   : lW_col   = mean - 2*sd. 'Left-Whisker' value (for boxplotting).
#   rW   : rW_col   = mean + 2*sd. 'Right-Whisker' value (for boxplotting).
#
#
# Note:
#   1)Place this script and get_FASTA_stats.py in the same working.direct
#   2)Scipy dependency
#   3)Numpy dependency
#   4)SeqIO from biopython 1.56 is used [note that the module is imported
#     when importing ger_FASTA_stats.py]
#
#
# Pau C. 01/04/2012
##############################################################################


import sys               # Collect arguments in command line
import os                # Interact with directory tree
import shutil            # Remove a directory
from scipy.stats import scoreatpercentile      # Use quartiles

# Use statistical functions and array creation and manipulation
from numpy import zeros, add, mean, median, subtract, std

# When we import get_FASTA_stats.py, a function [basic_stats()] is loaded,
from get_FASTA_stats import *
##############################################################################

# First we execute the function basic_stats() from the get_FASTA_stats.py file
# to get the three variables holding:
#    max_seq_len ->  [integer] length of longest sequence.
#    num_seq ->      [integer] total number of sequences.
#    nuc_len ->      [list]    length of each sequences.

(max_seq_len, num_seq, nuc_len) = basic_stats(fasta_file)

# Initialize vectors to hold bases count
A_count_partial = zeros(max_seq_len, int)
C_count_partial = zeros(max_seq_len, int)
G_count_partial = zeros(max_seq_len, int)
T_count_partial = zeros(max_seq_len, int)
N_count_partial = zeros(max_seq_len, int)

# Open the FASTA file. Rise an error if not possible.
try:
    handle = open(fasta_file, 'rU')
except IOError:
    print 'cannot open', fasta_file,
    '\nMake sure the file is in the working directory and you \
    have typed it correctly'

# Iterate over only the FASTA file with Bio.SeqIO.parse
# We center all sequences in its middle base, being
# len(my_seq)/2 the middle position.
# Note that odd values divided by 2 return just the integer part as in:
# 5/2 = 2 (2 comes from 2.5)

middle_longest = max_seq_len / 2  # Middle position of longest seq

for record in SeqIO.parse(handle, "fasta"):
        move_position = 0

        # len(record.seq) evaluates length of each sequence
        position = middle_longest - len(record.seq) / 2

        # Here we iterate over each base in the sequence and add 1
        for base in record.seq:
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
            # Here we capture bases that are not an A, C, G, T or N
            else:
                print "CAUTION!! Unexpected base <<" + base + ">> found\n"
            move_position += 1

# add() function from numpy is used to calculate total number of bases
# in each position
calc_matrix = add(A_count_partial, C_count_partial)
calc_matrix = add(calc_matrix, G_count_partial)
calc_matrix = add(calc_matrix, T_count_partial)
pos_total_bases = add(calc_matrix, N_count_partial)

handle.close()

# Open the QUAL file
# Note that 'your_file.fasta' is used to get the name of the QUAL file,
# therefore, the FASTA and QUAL file have to have the same name

file_name = fasta_file.split(".")[0]      # Get the name of the fasta file
qual_file = '.'.join([file_name, "qual"])  # Join two words by the separator
                                          # specified

# Check that is possible to open the QUAL file. Rise an error if not.
try:
    handle = open(qual_file, 'rU')
except IOError:
    print 'cannot open', qual_file,
    '\nMake sure the FASTA and QUAL files have the same name and\
    you have typed it correctly'
handle.close()

# Create a directory where temporary quality files will be placed.
# The strategy used is that each column in the aligned sequences,
# is going to be a separeted file which is going to be read afterwards
ini_wd = os.getcwd()
Qual_temp_dir = ini_wd + '/Qual_temp'
os.makedirs(Qual_temp_dir)


fasta_entry = 0
for record in SeqIO.parse(qual_file, "qual"):
    # Get the position where each sequence is centered
    current_len = nuc_len[fasta_entry]
    ini_centr_pos = middle_longest - current_len / 2

    fasta_entry += 1

    next_qual = 0
    for qual_val_pos in range(0, current_len):
        # Access each qual val one at a time
        qual_temp = record.letter_annotations['phred_quality'][next_qual]

        # Create the name of the file holding
        # Qual values for the given centerd column
        name = './Qual_temp/file_' + str(ini_centr_pos)

        handle = open(name, 'a')
        handle.write(str(qual_temp) + '\n')

        ini_centr_pos += 1
        next_qual += 1
handle.close()


# Calculate statistics over the temp quality files
# opening one at a time and reading the values
##################################################

min_col = zeros(max_seq_len, int)
max_col = zeros(max_seq_len, int)
sum_col = zeros(max_seq_len, int)
mean_col = zeros(max_seq_len, float)
Q1_col = zeros(max_seq_len, float)
med_col = zeros(max_seq_len, float)
Q3_col = zeros(max_seq_len, float)
IQR_col = zeros(max_seq_len, float)
lW_col = zeros(max_seq_len, float)
rW_col = zeros(max_seq_len, float)

# Iterate as many times as columns in the longest sequence
for number in range(0, max_seq_len):
    file_qual_name = './Qual_temp/file_' + str(number)

    handle = open(file_qual_name, 'r')
    temp_qual_list = []

    for line in handle:
        line_int = int(line)
        temp_qual_list.append(line_int)

    # Calculate minimum value and place it into min_col
    min_col[number] = min(temp_qual_list)
    # Calculate maximum value and place it into max_col
    max_col[number] = max(temp_qual_list)
    # Calculate the addition of all values and place it into sum_col
    sum_col[number] = sum(temp_qual_list)
    # Calculate mean value and place it into mean_col
    mean_tmp = round(mean(temp_qual_list), 2)
    mean_col[number] = mean_tmp
    # Calculate median value and place it into med_col
    med_col[number] = median(temp_qual_list)
    # Calculate first quartile (Q1) and place it into Q1_col
    Q1_col[number] = scoreatpercentile(temp_qual_list, 25)
    # Calculate third quartile (Q3) and place it into Q3_col
    Q3_col[number] = scoreatpercentile(temp_qual_list, 75)
    # Calculate the standard deviation
    std_col_tmp = std(temp_qual_list)
    # Calculate left-whisker and place it into lW_col
    # By convention, left-whisker = mean - 2* std
    lW_col[number] = round(mean_tmp - 2 * std_col_tmp, 2)
    # Calculate right whisker and place it into rW_col
    # By convention, left-whisker = mean + 2* std
    rW_col[number] = round(mean_tmp + 2 * std_col_tmp, 2)

# Calculate IQR and place it into IQR_col
IQR_col = subtract(Q3_col, Q1_col)

handle.close()

# Remove the temporary directory
shutil.rmtree(Qual_temp_dir)
##########################################################################
# The following section deals with the edition of a tabulated file containig
# the statistics we just calculated

# Where the data is going to be writen
os.chdir(ini_wd)
file_to_write = "tab_table_" + file_name + '.tsv'
handle_out = open(file_to_write, 'w')

# Start writing data
handle_out.write('column\tcount\tmin\tmax\tsum\tmean\tQ1\tmed\tQ3\tIQR\tlW\t\
rW\tA_Count\tC_Count\tG_Count\tT_Count\tN_Count\tMax_count\n')

iter = 1
for number in range(0, max_seq_len):
    handle_out.write(str(iter))     # column
    handle_out.write('\t')
    handle_out.write(str(pos_total_bases[number]))  # count
    handle_out.write('\t')
    handle_out.write(str(min_col[number]))
    handle_out.write('\t')
    handle_out.write(str(max_col[number]))
    handle_out.write('\t')
    handle_out.write(str(sum_col[number]))
    handle_out.write('\t')
    handle_out.write(str(mean_col[number]))
    handle_out.write('\t')
    handle_out.write(str(Q1_col[number]))
    handle_out.write('\t')
    handle_out.write(str(med_col[number]))
    handle_out.write('\t')
    handle_out.write(str(Q3_col[number]))
    handle_out.write('\t')
    handle_out.write(str(IQR_col[number]))
    handle_out.write('\t')
    handle_out.write(str(lW_col[number]))
    handle_out.write('\t')
    handle_out.write(str(rW_col[number]))
    handle_out.write('\t')
    handle_out.write(str(A_count_partial[number]))
    handle_out.write('\t')
    handle_out.write(str(C_count_partial[number]))
    handle_out.write('\t')
    handle_out.write(str(G_count_partial[number]))
    handle_out.write('\t')
    handle_out.write(str(T_count_partial[number]))
    handle_out.write('\t')
    handle_out.write(str(N_count_partial[number]))
    handle_out.write('\t')
    handle_out.write(str(num_seq))      # Max_count
    handle_out.write('\n')
    iter += 1
handle_out.close()
