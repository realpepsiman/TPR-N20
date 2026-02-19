# TPR-N20
Contains the script used to filter TPR-N20 nanopore sequencing data.
Takes a fastq sequencing file and a Levenshtein distance threshold score that determines how many changes can be present in the aptamer sequence. Set to 0 for only the perfect aptamer sequence.
So far, I have used 10 as the threshold which works fine.

The script creates a new fasta file with the suffix _filtered.fasta, which only contains the sequences with a subsequence "close enough" to the aptamer. The sequencing reads that are antisense to the functional RNA have been reverse complemeted, have the prefix "rev" in their name and are found below the sense sequence reads.

# Dependencies
The Levenshtein python package

# Usage
Download 'aptamer_filter.py' from this page.
In the command prompt, type:
`$python <filename.fastq> <threshold>`

# Example
`$python XX3SJJ_2_sample_2_PCR.fastq 10`
