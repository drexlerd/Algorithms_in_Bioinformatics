"""
Module name: Needleman Wunsch
Module author: dominik drexler <drexlerd@informatik.uni-freiburg.de>
"""

from prakt.nw import NeedlemanWunschBase
from prakt.fasta_parser.fasta_parser import parse_fasta, check_sequences_alphabet, SequenceType
from prakt.scoring_func_parser.scoring_func_parser import ScoringMatrix, MetricType
from prakt.basis_classes.cell import Cell
from prakt.util.util import compute_traceback
from enum import Enum
import argparse


class Case(Enum):
    ALIGNED = 0,
    SEQ1_GAPPED = 1,
    SEQ2_GAPPED = 2,


class Info(Enum):
    OK = 0,
    WRONG_ALPHABET = 1,


@NeedlemanWunschBase.register
class NeedlemanWunsch(NeedlemanWunschBase):
    """Document me!"""

    def base_case_init(self, d, seq1, seq2, scoring_matrix):
        """The initialization of the dynamic programming matrix d

        Args:
          d (list(list(Cell))): The DP matrix of correct shapes
          seq1 (str): The first sequence
          seq2 (str): The second sequence
          scoring_matrix (ScoringMatrix): The scoring matrix
        """

        # base case 
        for i in range(1, len(seq1) + 1):
            d[i][0].SetValue(i * scoring_matrix.cost_gap_open)
            d[i][0].AddPredecessor(d[i-1][0], Case.SEQ2_GAPPED)
        for j in range(1, len(seq2) + 1):
            d[0][j].SetValue(j * scoring_matrix.cost_gap_open)
            d[0][j].AddPredecessor(d[0][j-1], Case.SEQ1_GAPPED)


    def fill_matrix(self, d, seq1, seq2, scoring_matrix):
        """Fills out the dynamic programming matrix for the given seq1 and seq2
        with given scoring function and function_operation

        Args:
          d (list(list(Cell))): The preinitialized DP matrix
          seq1 (string): first sequence
          seq2 (string): second sequence
          scoring_matrix (ScoringMatrix): the scoring function object
        """
        for i in range(1, len(seq1) + 1):
            for j in range(1, len(seq2) + 1):
                # fill the matrix
                x_i = seq1[i-1]
                y_j = seq2[j-1]
                
                if x_i == "X" or y_j == "X":  # used in feng doolittle
                    sequence = [(d[i-1][j-1].value, Case.ALIGNED), 
                                (d[i][j-1].value, Case.SEQ1_GAPPED), 
                                (d[i-1][j].value, Case.SEQ2_GAPPED)]
                else:
                    score_aligning = scoring_matrix.score(x_i, y_j)
                    sequence = [(d[i-1][j-1].value + score_aligning, Case.ALIGNED), 
                                (d[i][j-1].value + scoring_matrix.cost_gap_open, Case.SEQ1_GAPPED), 
                                (d[i-1][j].value + scoring_matrix.cost_gap_open, Case.SEQ2_GAPPED)]
                result_sequence = scoring_matrix.function_operation(sequence)
                # store result of the recursion
                d[i][j].SetValue(result_sequence[0][0])  # set value
                # set predecessors
                for result_case in result_sequence:
                    if result_case[1] == Case.ALIGNED:  # aligning x_i with y_j
                        d[i][j].AddPredecessor(d[i-1][j-1], Case.ALIGNED)
                    elif result_case[1] == Case.SEQ1_GAPPED:  # gap in seq1
                        d[i][j].AddPredecessor(d[i][j-1], Case.SEQ1_GAPPED)
                    elif result_case[1] == Case.SEQ2_GAPPED:  # gap in seq2
                        d[i][j].AddPredecessor(d[i-1][j], Case.SEQ2_GAPPED)


    def traceback(self, d, seq1, seq2, complete_traceback=False):
        """Given the dynamic programming matrix and the two sequences
        this function computes all optimal alignments and returns them
        or just 1 if the complete_traceback flag is set.

        Args:
          d (list(list(Cell))): The finished DP matrix with values and predecessor pointers
          seq1 (str): The first sequence
          seq2 (str): The second sequence
          complete_traceback (bool): If True, returns all tracebacks, else 1

        Returns:
          list(list(str)): A list of optimal pairwise alignments
        """
        tracebacks = compute_traceback(d[-1][-1], all)
        alignments = []
        for traceback in tracebacks:
            alignment_seq1 = ""
            alignment_seq2 = ""
            # current positions in the sequences
            i = 0
            j = 0
            for node, case in reversed(traceback):
                if case == Case.ALIGNED:  # aligned
                    alignment_seq1 += seq1[i]
                    alignment_seq2 += seq2[j]
                    i += 1
                    j += 1
                elif case == Case.SEQ1_GAPPED:  # gap in seq1
                    alignment_seq1 += "_"
                    alignment_seq2 += seq2[j]
                    j += 1
                elif case == Case.SEQ2_GAPPED:  # gap in seq2
                    alignment_seq1 += seq1[i]
                    alignment_seq2 += "_"
                    i += 1
            alignments.append((alignment_seq1, alignment_seq2))
        if not complete_traceback:
            return [alignments[0]]
        return alignments


    def compute_optimal_alignments(self, seq1, seq2, scoring_matrix, complete_traceback):
        """Computes optimal pairwise alignments between seq1 and seq2 with given scoring matrix

        Args:
          seq1 (str): The first sequence
          seq2 (str): The second sequence
          scoring_matrix (ScoringMatrix): The scoring matrix
          complete_traceback (bool): If True, returns all tracebacks, else 1

        Returns:
          list(list(str)) : A 2D-array containing information about the pairwise optimal alignments
        """

        # dp matrix
        d = [[Cell(i, j) for j in range(len(seq2) + 1)] for i in range(len(seq1) + 1) ]

        # base cases
        self.base_case_init(d, seq1, seq2, scoring_matrix)

        # recursive case
        self.fill_matrix(d, seq1, seq2, scoring_matrix)

        #print("d")
        #for i in range(len(seq1) + 1):
        #    for j in range(len(seq2) + 1):
        #        print("%3d" % (d[i][j].value), end='')
        #    print()

        alignments = self.traceback(d, seq1, seq2, complete_traceback)

        score = d[-1][-1].value

        return score, alignments


    def run(self,
            seq1_fasta_fn,
            seq2_fasta_fn,
            subst_matrix_fn,
            is_distance_fn,
            cost_gap_open,
            complete_traceback):
            """
            Compute all optimal pairwise alignments between the sequences
            in the given fasta file seq1_fasta_fn and seq2_fasta_fn

            Args:
              seq1_fasta_fn (str): The relative path to a fasta file
              seq2_fasta_fn (str): The relative path to a fasta file
              subst_matrix_fn (str): The relative path to a scoring matrix file
              is_distance_fn (bool): If True, handle scoring matrix as distance measure, else similarity measure
              cost_gap_open (int): gap cost open
              complete_traceback (bool): If True, returns all tracebacks, else 1

            Returns:
              list(list(str)) : A 2D-array containing information about the pairwise optimal alignments
            """
            # scoring function
            scoring_matrix = ScoringMatrix(subst_matrix_fn, is_distance_fn, cost_gap_open)

            # sequences with their ids
            records_f1 = parse_fasta(seq1_fasta_fn)
            records_f2 = parse_fasta(seq2_fasta_fn)

            # check if the sequences are legal
            if not check_sequences_alphabet(records_f1, SequenceType.PROTEIN) \
              or not check_sequences_alphabet(records_f2, SequenceType.PROTEIN):
                return None, Info.WRONG_ALPHABET
            

            # init result array
            result = [[None for _ in records_f2] for _ in records_f1]

            for i in range(len(records_f1)):
                record1 = records_f1[i]
                for j in range(len(records_f2)):
                    record2 = records_f2[j]
                    seq1 = str(record1.seq)
                    seq2 = str(record2.seq)
                    score, alignments = self.compute_optimal_alignments(seq1, seq2, scoring_matrix, complete_traceback)
                    result[i][j] = (record1, record2, alignments, score) 

            return result, Info.OK


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Needleman Wunsch command line tool")
    parser.add_argument("seq1_fasta_fn", type=str)
    parser.add_argument("seq2_fasta_fn", type=str)
    parser.add_argument("subst_matrix_fn", type=str)
    parser.add_argument("cost_gap_open", type=int)
    parser.add_argument("--d", "--is_distance_fn", action='store_true')
    parser.add_argument("--c", "--complete_traceback", action='store_true')
    args = parser.parse_args()
    # run Needleman-Wunsch with some parameters
    nw = NeedlemanWunsch()

    result, info = nw.run(
        args.seq1_fasta_fn,
        args.seq2_fasta_fn,
        args.subst_matrix_fn,
        args.d,
        args.cost_gap_open,
        args.c)

    if info == Info.WRONG_ALPHABET:
        print("Wrong alphabet in sequences")
        exit(1)

    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print("Needleman Wunsch - Results")
    print("Scoring function: %s" % (args.subst_matrix_fn))
    print("Scoring type: %s" % ("Distance" if args.d else "Similarity"))
    print("Total amount of optimal aligmments: %d" % (len(result)))
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    seqs1_size = len(result)
    seqs2_size = len(result[0])
    for i in range(seqs1_size):
        for j in range(seqs2_size):
            print("Optimal pairwise alignments of sequences S_%d and S_%d" % (i+1,seqs1_size+j+1))
            score = result[i][j][3]
            print("Optimal Score: %.2f" % (score))
            for a in result[i][j][2]:
                print(a[0])
                print(a[1])
                print()
