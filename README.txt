PyEdFastX is a set of scripts that manipulate a FASTA and a QUAL file 
(it aims at Sanger technology, however, any FASTA and QUAL can be used)
to produce a table with some summary statistics ready to be read by 
FASTX-Toolkit (read more here:http://hannonlab.cshl.edu/fastx_toolkit/)

Pau Corral
########################################################################
Content of this README:

0)Files/scripts needed
1)Dependencies
2)Overview of what these scripts do
4)Usage example

########################################################################

0)Files/scripts needed
***********************
scripts:
    get_FASTA_stats.py
    table_FASTA_stats.py
files:
    your_file.fasta
    your_file.qual
(note that these .fasta and .qual files follow standard conventions on 
information layout and same file name)

1)Dependencies
***************
All these packages were downloaded from the ubuntu repository
Biopython 1.56 or above
Scipy 0.8.0 or above
Numpy 1.5.1 or above
Python (I used 2.6.5)

2)Overview of what these scripts do
************************************
Essentialy, PyEdFastX is a project (Python Editing for FastX) that contains two Python 
scripts that retrive some statistics from a given FASTA and QUAL files and organizes 
the output in a .tsv file (teb delimited file). One assumption that affects all the rest 
statistics is that all the sequences are aligned after fixing the middle position of all 
of them and spanning to both sides from this fixed position. Strictly speacking, this is not a
phylogenetic alignment, what we aim is to compare quality values in the extremes and 
in the middle positions, and to do such comparisons, the middle position of each sequence is aligned.

The statistics that the program reports are:

base positions          # counting according to longest sequence  
bases found in each position      # Counts of A, C, T and G
A base counts in each position
C base counts in each position
G base counts in each position
T base counts in each position
N base counts in each position
number of sequences in total in the file
Lowest quality score value found a given column
Highest quality score value found in a given column
Sum of quality score values for a given column
Mean quality score value for a given column
1st quartile quality score for a given column
Median quality score for a given column
3rd quartile quality score for a given column
Inter-Quartile range (Q3-Q1) for a given column
mean - 2*sd. 'Left-Whisker' value (for boxplotting)
mean + 2*sd. 'Right-Whisker' value (for boxplotting)

table_FASTA_stats.py is executed in the command line, passing the name of the 
fasta file to be analysed as an argument. A temporary folder (./Qual_temp) is 
created containg as many files (of the form 'file_#' where # is a number) 
as positions in the longest sequence. After the calculations, this directory is deleted. 

Finally, a .tsv file (of the form 'tab_table_' + your_fasta_file_name + '.tsv') 
is produced in the directpry where you first called the program.


3)Usage example
****************
Imagine you are going to analyse these files: orchid.fasta and orchid.qual containing
134000 sequences, the longest one having 900 nucleotides and the shortest 200. 
Working at bash shell with python interpreter installed, type:

$~> python table_FASTA_stats.py path_to_orchid.fasta

A correct execution would produce the file 'tab_table_orchid.tsv' in the 
directory where python was executed and the following output would be
printed into the screen:

Total number of sequences:  134000
Longest sequence: 	 900 nucl.
Shortest sequence: 	 200 nucl.

In that case, you have produced the desired .tsv file. 