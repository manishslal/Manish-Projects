import re
from matplotlib.pyplot import *
from Bio import pairwise2

def get_sequences_from_fasta(file_name):
    '''Get sequences from a fasta file. Returns a list of sequences.
    Example: get_sequences_from_file("seqs.fasta", "fasta")
    '''

    with open(file_name) as sequence_file:
        sequence_list = []
        sequence = ""
        for line in sequence_file:
            line = line.rstrip()  # strip end-of-line character from the end of the line
            if line.startswith(">"):
                # Get the first part of the header
                header = re.search("(>)(.+?)\s", line).group(2)
                if sequence:  # the presence of a header means the previous sequence record is finished. Store sequence into list now.
                    sequence_list.append(sequence)
                    sequence = ""
            else:
                sequence = sequence + line
        sequence_list.append(sequence)
        return sequence_list

def dot_plot(seq1, seq2, window=2):
    """make dot plot for two sequences.
    """

    # compare subsequences
    data = [[(seq1[i:i+window] != seq2[j:j+window])
                for j in range(len(seq1)-window)]
            for i in range(len(seq2)-window)]
    gray()
    xlabel("%s (length %i bp)" % ("rno-mir-150", len(seq1)))
    ylabel("%s (length %i bp)" % ("hsa-mir-150", len(seq2)))
    imshow(data)

def global_alignment(seq1, seq2, match, mismatch, gap_open, gap_extend):
    ''' implements Needleman-Wunsch global alignment
    '''
    alignments = pairwise2.align.globalms(seq1, seq2, match, mismatch, gap_open, gap_extend)
    for i in alignments:
        return format_alignment(*i)  # The splat operator unpacks a tuple into function arguments

def local_alignment(seq1, seq2, match, mismatch, gap_open, gap_extend):
    ''' implements Needleman-Wunsch global alignment
    '''
    alignments = pairwise2.align.localms(seq1, seq2, match, mismatch, gap_open, gap_extend)
    for i in alignments:
        return format_alignment(*i)  # The splat operator unpacks a tuple into function arguments

def format_alignment(align1, align2, score, begin, end):
      """Format the alignment prettily into a string."""
      s = []
      # s.append("%s\n" % align1)
      # s.append("%s%s\n" % (" " * begin, "|" * (end - begin)))
      # s.append("%s\n" % align2)
      s.append("Score=%g\n" % score)
      match = []
      for a, b in zip(align1, align2):
              if a == b:
                      match.append('|')
              else:
                      match.append('*')

      m="".join(match)
      s.append(align1+'\n')
      s.append(m+'\n')
      s.append(align2+'\n')

      return("".join(s))

def main():
    pass


if __name__ == '__main__':
    main()
    
# dot plot code modified from http://biopython.org/DIST/docs
