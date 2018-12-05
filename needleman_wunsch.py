from prakt.nw import NeedlemanWunschBase
from prakt.fasta_parser.fasta_parser import parse_fasta, check_sequences_alphabet, SequenceType
from prakt.scoring_func_parser.scoring_func_parser import ScoringMatrix, MetricType
from prakt.basis_classes.cell import Cell
from prakt.util.util import Min, Max, compute_traceback
from enum import Enum
import argparse


class Case(Enum):
    ALIGNED = 0,
    SEQ1_GAPPED = 1,
    SEQ2_GAPPED = 2,


@NeedlemanWunschBase.register
class NeedlemanWunsch(NeedlemanWunschBase):
    """Document me!"""

    def base_case_init(self, d, seq1, seq2, cost_gap_open):
        """The initialization of the dynamic programming matrix d

        Args:
          len_seq1 (int): The size of the first sequence
          len_seq2 (int): The size of the second sequence
        """

        # base case 
        for i in range(1, len(seq1) + 1):
            d[i][0].SetValue(i * cost_gap_open)
            d[i][0].AddPredecessor(d[i-1][0], Case.SEQ2_GAPPED)
        for j in range(1, len(seq2) + 1):
            d[0][j].SetValue(j * cost_gap_open)
            d[0][j].AddPredecessor(d[0][j-1], Case.SEQ1_GAPPED)


    def fill_matrix(self, d, seq1, seq2, scoring_matrix, function_operation, cost_gap_open):
        """Fills out the dynamic programming matrix for the given seq1 and seq2
        with given scoring function and function_operation

        Args:
          seq1 (string): first sequence
          seq2 (string): second sequence
          scoring_matrix (ScoringMatrix): the scoring function object
          function_operation (f(list(tuples))) : function returning all tuples with min/max tuple.item1
        """
        for i in range(1, len(seq1) + 1):
            for j in range(1, len(seq2) + 1):
                # fill the matrix
                x_i = seq1[i-1]
                y_j = seq2[j-1]
                score_aligning = scoring_matrix.score(x_i, y_j)

                sequence = [(d[i-1][j-1].value + score_aligning, Case.ALIGNED), 
                            (d[i][j-1].value + cost_gap_open, Case.SEQ1_GAPPED), 
                            (d[i-1][j].value + cost_gap_open, Case.SEQ2_GAPPED)]
                result_sequence = function_operation(sequence)
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
        """
        tracebacks = compute_traceback(d[-1][-1], all)
        alignments = []
        for traceback in tracebacks:
            alignment_seq1 = ""
            alignment_seq2 = ""
            # current positions in the sequences
            i = 0
            j = 0
            for case in reversed(traceback):
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


    def compute_optimal_alignments(self, seq1, seq2, scoring_matrix, cost_gap_open, complete_traceback):
        # scoring scheme (min / max)
        function_operation = None

        if scoring_matrix.metric_type == MetricType.DISTANCE:
            function_operation = Min
            cost_gap_open = abs(cost_gap_open)
        elif scoring_matrix.metric_type == MetricType.SIMILARITY:
            function_operation = Max
            cost_gap_open = - abs(cost_gap_open)

        # dp matrix
        d = [[Cell(i, j) for j in range(len(seq2) + 1)] for i in range(len(seq1) + 1) ]

        # base cases
        self.base_case_init(d, seq1, seq2, cost_gap_open)

        # recursive case
        self.fill_matrix(d, seq1, seq2, scoring_matrix, function_operation, cost_gap_open)

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

            Returns:
              List(List) : A 2D-array containing information about the pairwise optimal alignments
            """
            # scoring function
            scoring_matrix = ScoringMatrix(subst_matrix_fn, is_distance_fn)

            # sequences with their ids
            records_f1 = parse_fasta(seq1_fasta_fn)
            records_f2 = parse_fasta(seq2_fasta_fn)

            # check if the sequences are legal
            check_sequences_alphabet(records_f1, SequenceType.PROTEIN)
            check_sequences_alphabet(records_f2, SequenceType.PROTEIN)

            # init result array
            result = [[None for _ in records_f2] for _ in records_f1]

            for i in range(len(records_f1)):
                record1 = records_f1[i]
                for j in range(len(records_f2)):
                    record2 = records_f2[j]
                    seq1 = str(record1.seq)
                    seq2 = str(record2.seq)
                    score, alignments = self.compute_optimal_alignments(seq1, seq2, scoring_matrix, cost_gap_open, complete_traceback)
                    result[i][j] = (record1, record2, alignments, score) 

            return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Needleman Wunsch command line tool")
    parser.add_argument("seq1_fasta_fn", type=str)
    parser.add_argument("seq2_fasta_fn", type=str)
    parser.add_argument("subst_matrix_fn", type=str)
    parser.add_argument("is_distance_fn", type=bool)
    parser.add_argument("cost_gap_open", type=int)
    parser.add_argument("complete_traceback", type=bool)
    args = parser.parse_args()
    # run Needleman-Wunsch with some parameters
    nw = NeedlemanWunsch()
    print(nw.run(
        args.seq1_fasta_fn,
        args.seq2_fasta_fn,
        args.subst_matrix_fn,
        args.is_distance_fn,
        args.cost_gap_open,
        args.complete_traceback))
