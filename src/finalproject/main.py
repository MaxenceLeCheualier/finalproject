import argparse 
import sys 
from finalproject.fasta_reader import FastaReader
from finalproject.needleman_wunsch import NeedlemanWunsch 

def main():
    parser = argparse.ArgumentParser(description= "Sequence alignment using Needleman-Wunsch algorithm")
    
    parser.add_argument("-i", "--input", required = True, help = "Fasta file input")
    parser.add_argument("-o", "--output", default = "-", help = "Output file")
    parser.add_argument("-m", "--match-score", type = int, default = "1", help = "Score for a match")
    parser.add_argument("-x", "--mismatch-score", type = int, default = -1, help = "Score for a mismatch")
    parser.add_argument("-j", "--indel-score", type = int, default = -2, help = "Score for instertion / deletion")
    parser.add_argument("-v", "--verbose", action = "count", default = 0, help = "Verbosity level") 

    args = parser.parse_args()

    reader = FastaReader(args.input)
    sequences = reader.read()

    NW = NeedlemanWunsch(sequences, 
                                 match_score= args.match_score, 
                                 mismatch_penalty = args.mismatch_score, 
                                 gap_penalty = args.indel_score
                                 )
    
    alignement = NW.align(sequences)

    if args.output == "-":
        out = sys.stdout
    else:
        out = open(args.output, "w")

    for name, aligned in alignement.items():

        seq1 = aligned.original
        seq2 = aligned.variant

        print(f"{name}:", file=out)
        print(seq1, file=out)
        print(seq2, file=out)
        print("", file=out)

    if args.output != "-":
        out.close()


if __name__ == "__main__":
    main()

    