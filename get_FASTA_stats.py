#
# This scrip creates a funtion [basic_stats()] that calculates:
#    num_seq ->      [integer] how many sequences are in the fasta file.
#    max_seq_len ->  [integer] length of longest sequence.
#    min_seq_len ->  [integer] length of the shortest sequence.
#    nuc_len ->      [list]    length of each of the sequences.
#
# The function returns these three values:
# max_seq_len, num_seq, nuc_len
#
# Command to execute in the command line:
#    $~> python -i get_FASTA_stats file_name.fasta
# and when python interactive mode is entered, type:
#    >>>basic_stats(fasta_file)
# This will print into STDOUT and return  into the workspace:
# max_seq_len, num_seq, nuc_len
#
# SeqIO from biopython 1.56 (Ubuntu repository) is used to parse the FASTA
#
# Pau C.     1/04/2012
################################################################

from Bio import SeqIO                   # Use the parser
import sys                              # Collect arguments in command line

################################################################

fasta_file = sys.argv[1]  # Where, your_file.fasta = sys.argv[1]
                          # ex.: $~> python get_FASTA_stats.py your_file.fasta


def basic_stats(file):
    """Return two integers with total number of sequences and maximum length,
    and a list with the length of each sequence
    """
    nuc_len = []           # Sequence length - list
    max_seq_len = 0        # Longest sequence - integer
    num_seq = 0            # Total number of sequences - integer
    min_seq_len = 500      # Shortest sequence - integer
                           # set to 500 as it is the average length
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
        # Here the shortest
        elif nuc_length < min_seq_len:
            min_seq_len = nuc_length
    handle.close()

    print "Total number of sequences: ", num_seq
    print "Longest sequence: \t %i nucl." % max_seq_len
    print "Shortest sequence: \t %i nucl." % min_seq_len

    return  max_seq_len, num_seq, nuc_len
