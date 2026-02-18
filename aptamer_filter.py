# Takes a FASTQ file containing nanopore sequencing data from TPR-N20 directed evolution
# creates a FASTA file containing only the sequences that contain a sequence similar
# to the strep aptamer based on the threshold given
import sys
import Levenshtein

apt_seq = "GACCGACCAGAATCATGCAAGTGCGTAAGATAGTCGCGGGCCGGG"
apt_seq_rev = "CCCGGCCCGCGACTATCTTACGCACTTGCATGATTCTGGTCGGTC"
ref_seq = "GAAATTAATACGACTCACTATAGGATCTTCTCGATCTAACAAAAAAGACAAATCTGCCACAAAGCTTGAGAGCATCTTCGGATGCAGAGGCGGCAGCCTTCGGTGGCGCGATAGCGCCAACGTTCTCAACTATGACACGCAAAACGCGTGCTCCGTTGAATGGAGTTTATCATGNNNNNNNNNNNNNNNNNNNNAACAAACAAACAGACCGACCAGAATCATGCAAGTGCGTAAGATAGTCGCGGGCCGGGGGAATCGATGCCGAA"

def reverse_complement(seq):
    # takes a DNA sequence, returns the reverse complemented sequence in the 5'-3' direction
    complement = []
    reverse_dict = {'a': 't', 'A': 'T', 't': 'a', 'T': 'A', 'g': 'c', 'G': 'C', 'c': 'g', 'C': 'G', 'n': 'n', 'N': 'N', '\n': '\n'}
    for base in seq[::-1]:
        complement.append(reverse_dict[base])
    return ''.join(complement)


def convert_to_fasta(filename):
    # reads FASTQ file and converts it to FASTA format
    seqs = []
    f = open(filename, "r")
    f_enum = list(enumerate(f))
    # discard lines containing sequencing quality data (add only the first 2 lines from each sequence entry)
    # adds the currently indexed and next line in the file to the seqs list, then skips three indexes and repeats the process
    # also replaces @ with > to follow FASTA format
    line_n = 1
    for i, line in f_enum:
        if line_n == 1:
            seqs.append(f_enum[i][1].replace("@", ">")) 
            seqs.append(f_enum[i+1][1])
            line_n +=1
        elif line_n == 4:
            line_n = 1
        else:
            line_n +=1
    f.close()
    return seqs

def filter_levenshtein(seqs, threshold):
    # removes sequences without subsequence similar to aptamer (based on threshold)
    # reverse complements sequences if they contain a sequence similar to reverse complemented aptamer
    filtered_seqs = []
    filtered_seqs_rev = []
    for i in range(len(seqs)):
        high_score = 1000
        high_score_rev = 1000 # arbitrarily high number, score should always be lower (which is better)
        for j in range(len(seqs[i])-50):
            substring = seqs[i][j:j+50]
            score = Levenshtein.distance(substring, apt_seq)
            if  score < high_score:
                high_score = score
            score_rev = Levenshtein.distance(substring, apt_seq_rev)
            if score_rev < high_score_rev:
                high_score_rev = score_rev
        if high_score < threshold:
            filtered_seqs.append(seqs[i-1])
            filtered_seqs.append(seqs[i])
        if high_score_rev < threshold:
            filtered_seqs_rev.append(">rev-" + seqs[i-1][1:]) # name of sequence is modified to denote that this sequence has been reversed
            filtered_seqs_rev.append(reverse_complement(seqs[i]))
    filtered_seqs += filtered_seqs_rev # reversed sequences will appear below non-reversed sequences
    return filtered_seqs
        

def create_fasta(filename, seqs):
    # creates a FASTA file containing sequences using the FASTQ filename
    f = open(filename[:-6] + "_filtered.fasta", 'w')
    for line in seqs:
        f.write(line + '\n')
    f.close

def filter_seqs(filename, threshold):
    # does the thing
    seqs = convert_to_fasta(filename)
    filtered_seqs = filter_levenshtein(seqs, threshold)
    create_fasta(filename, filtered_seqs)


# command line argument stuff
if __name__ == "__main__":
  filename = sys.argv[1]
  threshold = int(sys.argv[2])
  filter_seqs(filename, threshold)
