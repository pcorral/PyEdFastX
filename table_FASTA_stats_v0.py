#
# Produce a table with statistics from a FASTA & QUAL file
# still incomplete...
#
# Essentially, all the preprocessint deals with the fact that what we are trying to do is
# alaign all sequences by its middle base, so that, both extrems in each sequence can be compared
#
#
# Note:
#	1)place get_FASTA_stats.py and table_FASTA_stats.py in the same working directory.
#	2)execute get_FASTA_stats.py with basic_stats() funtion to know
#	  the length of the longest sequence
#
#
# Pau C. (initiated 01/04/2012)
##################################################


import sys									#collect arguments in command line

from scipy.stats import scoreatpercentile					#use quartiles
from Bio.SeqIO.QualityIO import PairedFastaQualIterator 	#use FASTA & QUAL parsers
from numpy import *											#use statistic functions as sum, min, max, mean
															#use array creation and manipulation
#from Bio import SeqIO 						#this is loaded when get_FASTA_stats.py is executed											


#When we import this file, a function [basic_stats()] is executed, and
#two values are loaded into the workspace:
#	->length of the longest sequence
#	->number of sequences

from get_FASTA_stats import *		
###################################################################

#first we execute the function basic_stats() from the get_FASTA_stats.py file
#to get the two variables holding
#	longest seq in the FASTA -> longest_out
#	number of seq's in the FASTA -> total_seq_out

(longest_out, total_seq_out) = basic_stats(file)

#Initialize vectors to hold bases count
A_count_partial = zeros(longest_out,int)
C_count_partial = zeros(longest_out,int)
G_count_partial = zeros(longest_out,int)
T_count_partial = zeros(longest_out,int)
N_count_partial = zeros(longest_out,int)

#Open the FASTA file. Rise an error if not possible.
try:
	handle = open(file, 'rU')
except IOError:
	print 'cannot open', file, 
	'\nMake sure the file is in the working directory and you have typed it correctly'

#Iterate over only the FASTA file with the Biopython parser in SeqIO
#We center all sequences in its middle base, being len(my_seq)/2 the middle position

for record in SeqIO.parse(handle, "fasta"):
	
	move_position = -1  				# Increment +1 in each iteration
										# set to -1 to start the iteration with value 0
										# (0 points to first position)
	middle_longest = longest_out/2					# middle position for the longest sequence

	position = middle_longest - len(record.seq)/2	# len(record.seq) evaluates length of each sequence
													# note that odd values divided by 2 return just the integer part as in:
													# 5/2 = 2 (2 comes from 2.5)
	
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
	
#add() function from numpy is used to calculate total number of bases in each position
calc_matrix = add(A_count_partial, C_count_partial)
calc_matrix = add(calc_matrix, G_count_partial)
calc_matrix = add(calc_matrix, T_count_partial)
per_position_total_bases = add(calc_matrix, N_count_partial)

print per_position_total_bases 

# What I have so far is: (which has to be writen in the table [tabulated file])
# Nomenclature -->
#   in_table_name : in_my_code_name = explanation of what it is/ what it holds
#
#	column	: extremes 					= range of numbers in the longest sequence from the first position till the last
#	count	: per_position_total_bases 	= count of all bases found in each position
#	A_Count	: A_count_partial 			= A base counts in each position
#	C_Count : C_count_partial			= C base counts in each position
#	G_Count : G_count_partial 			= G base counts in each position
#	T_Count : T_count_partial 			= T base counts in each position
#	N_Count : N_count_partial 			= N base counts in each position
#	Max_count : total_seq_out			= number of sequences
#
#still pending ....

#They only work with the QUAL values:
#min     = Lowest quality score value found in this column.
#max     = Highest quality score value found in this column.
#sum     = Sum of quality score values for this column.
#mean    = Mean quality score value for this column.
#Q1	= 1st quartile quality score.
#med	= Median quality score.
#Q3	= 3rd quartile quality score.
#IQR	= Inter-Quartile range (Q3-Q1).
#lW	= -2*sd. 'Left-Whisker' value (for boxplotting).
#rW	= 2*sd. 'Right-Whisker' value (for boxplotting).
