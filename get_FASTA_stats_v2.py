#New line added
#
# Parse a FASTA file with the defined function basic_stats() to abtain
#    num_seq - how many sequences are there. Value is saved as an integer
#    max_seq_len - maximum sequence length. Value is saved as an integer
# The function returns these two values
#
# command to execute:
#    $~> python -i get_FASTA_stats file_name_of_FASTA.fasta
#
# SeqIO from biopython 1.56 (Ubuntu repository) is used to parse the FASTA
#
#
# Pau C.     1/04/2012
################################################################

from Bio import SeqIO                     # use FASTA parser
import sys                                # collect arguments in command line


################################################################

fasta_file = sys.argv[1]      # where, your_file.fasta = sys.argv[1]
                          # ex.: $~> python get_FASTA_stats.py your_file.fasta


def basic_stats(file):
    """Return two integers with total number of sequences and maximum length, 
    and a list with the length of each sequence
    """
    nuc_len = []           # sequence length - list
    max_seq_len = 0        # longest sequence - integer
    num_seq = 0            # total number of sequences - integer

    # Try to open the FASTA file. If not, an error is raised.
    try:
        handle = open(file, 'rU')
    except IOError:
        print 'cannot open', file,
        "\nMake sure the file is in the working directory and \
         you have typed it correctly"

    # Read the fasta file, one sequence at a time
    for record in SeqIO.parse(handle, "fasta"):
        num_seq += 1
        nuc_length = len(record.seq)        # record.seq evaluates all nucleot.
        nuc_len.append(nuc_length)
        # Here we evaluate the longest sequence
        if nuc_length > max_seq_len:
            max_seq_len = nuc_length
    handle.close()

    print "Total number of sequences: \t", num_seq
    print "Longest sequence: \t\t", max_seq_len

    return  max_seq_len, num_seq, nuc_len
