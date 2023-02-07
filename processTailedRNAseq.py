import gzip
import sys
from Bio import SeqIO
from Bio import Align
from Bio.SeqIO.QualityIO import FastqGeneralIterator

# Arguments
# print("Arguments set!")
file1 = sys.argv[1]
file2 = sys.argv[2]
adaptor_root = sys.argv[3]
adaptor_five = sys.argv[4]
adaptor_three = sys.argv[5]
out_five = sys.argv[6]
out_three = sys.argv[7]
out_else = sys.argv[8]

output_five_file_1 = out_five + "_1.fastq"
output_five_file_2 = out_five + "_2.fastq"
output_three_file_1 = out_three + "_1.fastq"
output_three_file_2 = out_three + "_2.fastq"
output_else_file_1 = out_else + "_1.fastq"
output_else_file_2 = out_else + "_2.fastq"

output_five_1 = open(output_five_file_1, 'wt', 1000000)
output_five_2 = open(output_five_file_2, 'wt', 1000000)
output_three_1 = open(output_three_file_1, 'wt', 1000000)
output_three_2 = open(output_three_file_2, 'wt', 1000000)
output_else_1 = open(output_else_file_1, 'wt', 1000000)
output_else_2 = open(output_else_file_2, 'wt', 1000000)
# print("Write file handles open set!")

#Setup aligner with match 1, mismatch 0, local
aligner = Align.PairwiseAligner()
aligner.mode = 'local'
aligner.open_gap_score = 0


with open(file1, "rt") as handle1:
    with open(file2, "rt") as handle2:
        #print("Read files opened succesfully!")
        for (title1, seq1, qual1), (title2, seq2, qual2) in zip(FastqGeneralIterator(handle1), FastqGeneralIterator(handle2)):
            root_alignment = aligner.align(seq2[0:30], adaptor_root)
            if root_alignment.score > 20:
                five_alignment = aligner.align(seq2[0:30], adaptor_five).score
                three_alignment = aligner.align(seq2[0:30], adaptor_three).score
                # if a 5 prime alignment is likelier, then output it there, else output to 3 prime file
                if five_alignment > three_alignment:
                    output_five_1.write("@%s\n%s\n+\n%s\n" % (title1, seq1, qual1))
                    #These next lines clip all initial G's
                    #j = 30 + len(seq2[30:151]) - len(seq2[30:151].lstrip('G'))
                    #output_five_2.write("@%s\n%s\n+\n%s\n" % (title2, seq2[j:151], qual2[j:151]))
                    #This line doesnt, only clips the ones that are within the adaptor in the first 30 bases
                    output_five_2.write("@%s\n%s\n+\n%s\n" % (title2, seq2[30:151], qual2[30:151]))
                elif three_alignment > five_alignment:
                    output_three_1.write("@%s\n%s\n+\n%s\n" % (title1, seq1, qual1))
                    output_three_2.write("@%s\n%s\n+\n%s\n" % (title2, seq2[30:151], qual2[30:151]))
            else:
                output_else_1.write("@%s\n%s\n+\n%s\n" % (title1, seq1, qual1))
                output_else_2.write("@%s\n%s\n+\n%s\n" % (title2, seq2, qual2))


print("File processed succesfully!")

output_five_1.close()
output_five_2.close()
output_three_1.close()
output_three_2.close()
output_else_1.close()
output_else_2.close()
