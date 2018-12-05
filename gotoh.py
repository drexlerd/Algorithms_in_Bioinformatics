from prakt.gt import GotohBase
from prakt.fasta_parser.fasta_parser import parse_fasta
from prakt.scoring_func_parser.scoring_func_parser import ScoringMatrix, MetricType
from prakt.basis_classes.cell import Cell
from prakt.util.util import Min, Max, compute_traceback
from enum import Enum


class Case(Enum):
    ALIGN_D = 0
    JUMP_P = 1
    JUMP_Q = 2
    GAP_SEQ1_D = 3  # opening gap in seq1
    GAP_SEQ1_Q = 4  # extending gap in seq1
    GAP_SEQ2_D = 5  # opening gap in seq2
    GAP_SEQ2_P = 6  # extending gap in seq2


@GotohBase.register
class Gotoh(GotohBase):
    """Document me!"""


    def base_case_init(self, d, p, q, 
      seq1, seq2, 
      affine_cost_gap_open, 
      affine_cost_gap_extend, 
      extreme_value):        
        """The inplace initialization of the dynamic programming matrix d, p, q

        Args:
          d (matrix): DP matrix ending with alignment
          p (matrix): DP matrix ending with gap in seq2
          q (matrix) DP matrix ending with gap in seq1
          seq1 (string): The first sequence
          seq2 (string): The second sequence
          affine_cost_gap_open (int): Cost for opening a new gap
          affine_cost_gap_extend (int): Cost for extending a gap
          extreme_value (float): Maximum value used in recursion initialization, either infty or -infty
        """

        # base case 
        for i in range(1, len(seq1) + 1):
            if i == 1:
                d[i][0].SetValue(affine_cost_gap_open)
            else:
                d[i][0].SetValue(affine_cost_gap_open + (i-1) * affine_cost_gap_extend)
            d[i][0].AddPredecessor(d[i-1][0], Case.GAP_SEQ2_D)
            q[i][0].SetValue(extreme_value)
        for j in range(1, len(seq2) + 1):
            if j == 1:
                d[0][j].SetValue(affine_cost_gap_open)
            else:
                d[0][j].SetValue(affine_cost_gap_open + (j-1) * affine_cost_gap_extend)
            d[0][j].AddPredecessor(d[0][j-1], Case.GAP_SEQ1_D)
            p[0][j].SetValue(extreme_value)


    def fill_matrix(self, d, p, q, 
      seq1, seq2, 
      scoring_matrix, 
      function_operation, 
      affine_cost_gap_open, 
      affine_cost_gap_extend):
        """The inplace filling of the dynamic programming matrix d, p, q

        Args:
          d (matrix): DP matrix ending with alignment
          p (matrix): DP matrix ending with gap in seq2
          q (matrix) DP matrix ending with gap in seq1
          seq1 (string): The first sequence
          seq2 (string): The second sequence
          scoring_matrix (ScoringMatrix): The substitution matrix with its type (Similarity/Distance)
          function_operation: Min-Maximization depending on Distance/Similarity
          affine_cost_gap_open (int): Cost for opening a new gap
          affine_cost_gap_extend (int): Cost for extending a gap
        """
        for i in range(1, len(seq1) + 1):
            for j in range(1, len(seq2) + 1):
                x_i = seq1[i-1]
                y_j = seq2[j-1]
                score_aligning = scoring_matrix.score(x_i, y_j)

                # compute p
                sequence_p = [(d[i-1][j].value + affine_cost_gap_open, Case.GAP_SEQ2_D), 
                                (p[i-1][j].value + affine_cost_gap_extend, Case.GAP_SEQ2_P)]
                result_sequence_p = function_operation(sequence_p)
                p[i][j].SetValue(result_sequence_p[0][0])  # set value
                # compute q
                sequence_q = [(d[i][j-1].value + affine_cost_gap_open, Case.GAP_SEQ1_D), 
                                (q[i][j-1].value + affine_cost_gap_extend, Case.GAP_SEQ1_Q)]
                result_sequence_q = function_operation(sequence_q)
                q[i][j].SetValue(result_sequence_q[0][0])  # set value
                # compute d
                sequence_d = [(d[i-1][j-1].value + score_aligning, Case.ALIGN_D), 
                                (p[i][j].value, Case.JUMP_P), 
                                (q[i][j].value, Case.JUMP_Q)]
                result_sequence_d = function_operation(sequence_d)
                d[i][j].SetValue(result_sequence_d[0][0])  # set value

                # set predecessors
                for result_case in result_sequence_d:
                    if result_case[1] == Case.ALIGN_D:
                        d[i][j].AddPredecessor(d[i-1][j-1], Case.ALIGN_D)
                    elif result_case[1] == Case.JUMP_P:
                        d[i][j].AddPredecessor(p[i][j], Case.JUMP_P)
                    elif result_case[1] == Case.JUMP_Q:
                        d[i][j].AddPredecessor(q[i][j], Case.JUMP_Q)
                # set predecessors of p
                for result_case in result_sequence_p:
                    if result_case[1] == Case.GAP_SEQ2_D:
                        p[i][j].AddPredecessor(d[i-1][j], Case.GAP_SEQ2_D)
                    elif result_case[1] == Case.GAP_SEQ2_P:
                        p[i][j].AddPredecessor(p[i-1][j], Case.GAP_SEQ2_P)
                # set predessors of q
                for result_case in result_sequence_q:
                    if result_case[1] == Case.GAP_SEQ1_D:
                        q[i][j].AddPredecessor(d[i][j-1], Case.GAP_SEQ1_D)
                    elif result_case[1] == Case.GAP_SEQ1_Q:
                        q[i][j].AddPredecessor(q[i][j-1], Case.GAP_SEQ1_Q)
                

    def traceback(self, d, p, q, seq1, seq2, complete_traceback=True):
        """Compute traceback starting from bottom right cell in d

        Args:
          d (matrix): DP matrix ending with alignment (completely filled)
          p (matrix): DP matrix ending with gap in seq2 (completely filled)
          q (matrix) DP matrix ending with gap in seq1 (completely filled)
          seq1 (string): The first sequence
          seq2 (string): The second 
          complete_traceback (bool): if True, then all optimal alignments are returned, else just 1
        """
        tracebacks = compute_traceback(d[-1][-1], complete_traceback)
        alignments = []
        for traceback in tracebacks:
            alignment_seq1 = ""
            alignment_seq2 = ""
            # current positions in the sequences
            i = 0
            j = 0
            for case in reversed(traceback):  # dfs started at bottom right cell => need to reverse it
                #print(case)
                if case in [Case.ALIGN_D]:  # aligned
                    alignment_seq1 += seq1[i]
                    alignment_seq2 += seq2[j]
                    i += 1
                    j += 1
                elif case in [Case.GAP_SEQ2_D, Case.GAP_SEQ2_P]:
                    alignment_seq1 += seq1[i]
                    alignment_seq2 += "_"
                    i += 1
                elif case in [Case.GAP_SEQ1_D, Case.GAP_SEQ1_Q]:
                    alignment_seq1 += "_"
                    alignment_seq2 += seq2[j]
                    j += 1
            alignments.append((alignment_seq1, alignment_seq2))
            #print()
        if not complete_traceback:
            return [alignments[0]]
        return alignments


    def compute_optimal_alignments(self, seq1, seq2, 
      subst_matrix_fn, 
      is_distance_fn, 
      affine_cost_gap_open, 
      affine_cost_gap_extend, 
      complete_traceback):
        """Computes optimal alignment(s) between two given sequences.
        """
        # scoring function
        scoring_matrix = ScoringMatrix(subst_matrix_fn, is_distance_fn)

        # scoring scheme (min / max)
        function_operation = None
        extreme_value = None  # used in initialization of p and q

        if scoring_matrix.metric_type == MetricType.DISTANCE:
            function_operation = Min
            extreme_value = float("inf")
            affine_cost_gap_open = abs(affine_cost_gap_open)
            affine_cost_gap_extend = abs(affine_cost_gap_extend)
        elif scoring_matrix.metric_type == MetricType.SIMILARITY:
            function_operation = Max
            extreme_value = - float("inf")
            affine_cost_gap_open = - abs(affine_cost_gap_open)
            affine_cost_gap_extend = - abs(affine_cost_gap_extend)
        
        # ends with alignment
        d = [[Cell(i, j) for j in range(len(seq2) + 1)] for i in range(len(seq1) + 1) ]
        # ends with gap in seq2
        p = [[Cell(i, j) for j in range(len(seq2) + 1)] for i in range(len(seq1) + 1) ]
        # ends with gap in seq1
        q = [[Cell(i, j) for j in range(len(seq2) + 1)] for i in range(len(seq1) + 1) ]

        # base cases
        self.base_case_init(d, p, q, seq1, seq2, affine_cost_gap_open, affine_cost_gap_extend, extreme_value)

        # recursive case
        self.fill_matrix(d, p, q, seq1, seq2, scoring_matrix, function_operation, affine_cost_gap_open, affine_cost_gap_extend)

        #print("d")
        #for i in range(len(seq1) + 1):
        #    for j in range(len(seq2) + 1):
        #        print("%3d" % (d[i][j].value), end='')
        #    print()

        #print("p")
        #for i in range(len(seq1) + 1):
        #    for j in range(len(seq2) + 1):
        #        print("%3.0f" % (p[i][j].value), end='')
        #    print()

        #print("q")
        #for i in range(len(seq1) + 1):
        #    for j in range(len(seq2) + 1):
        #        print("%3.0f" % (q[i][j].value), end='')
        #    print()


        alignments = self.traceback(d, p, q, seq1, seq2, complete_traceback)

        score = d[-1][-1].value

        return score, alignments


    def run(self,
            seq1_fasta_fn,
            seq2_fasta_fn,
            subst_matrix_fn,
            is_distance_fn,
            affine_cost_gap_open,
            affine_cost_gap_extend,
            complete_traceback):
            # sequences with their ids
            records_f1 = parse_fasta(seq1_fasta_fn)
            records_f2 = parse_fasta(seq2_fasta_fn)

            # init result array
            result = [[None for _ in records_f2] for _ in records_f1]

            for i in range(len(records_f1)):
                record1 = records_f1[i]
                for j in range(len(records_f2)):
                    record2 = records_f2[j]
                    seq1 = str(record1.seq)
                    seq2 = str(record2.seq)
                    score, alignments = self.compute_optimal_alignments(seq1, seq2, subst_matrix_fn, is_distance_fn, affine_cost_gap_open, affine_cost_gap_extend, complete_traceback)
                    result[i][j] = (record1, record2, alignments, score)
            return result
