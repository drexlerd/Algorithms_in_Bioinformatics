from prakt.nussinov import NussinovBase
from prakt.fasta_parser.fasta_parser import parse_fasta, check_sequences_alphabet, SequenceType
from prakt.basis_classes.cell import Cell
from enum import Enum
from prakt.util.util import Max, compute_traceback

class Case(Enum):
    UNPAIRED = 0,
    PAIRED_SUCCESS = 1,
    PAIRED_FAIL = 2


@NussinovBase.register
class Nussinov(NussinovBase):
    def is_base_pair(self, c1, c2):
        if c1 == "A" and c2 == "U" or c1 == "U" and c2 == "A" \
          or c1 == "G" and c2 == "C" or c1 == "C" and c2 == "G":
            return True
        return False


    def fill_matrix(self, d, sequence, loop_length=1):
        for j_start in range(2, len(sequence)):
            i = 0
            j = j_start
            while True:
                # i and j are the indices iterated in diagonal
                # TODO
                sequence = [(d[i][j-1].value, d[i][j-1], Case.UNPAIRED)]
                
                for k in range(i, j - loop_length):
                    if self.is_base_pair(sequence[k],sequence[j]):
                        sequence.append((d[i][k-1].value + d[k+1][j-1].value + 1, d[k][j], Case.PAIRED_SUCCESS))
                    else:
                        sequence.append((0, d[k,j], Case.PAIRED_FAIL))

                result_sequence = Max(sequence)

                # store result of the recursion
                d[i][j].SetValue(result_sequence[0][0])  # set value
                # set predecessors
                for result_case in result_sequence:
                    d[i][j].AddPredecessor(result_case[1], result_case[2])

                # END TODO
                i += 1
                j += 1
                if i >= len(sequence) or j >= len(sequence):
                    break


    def compute_optimal_structure(self, sequence, complete_traceback):
        d = [[Cell(i, j) for j in range(len(sequence) + 1)] for i in range(len(sequence))]

        self.fill_matrix(d, sequence, loop_length=1)

        tracebacks = compute_traceback(d[0][-1], all=complete_traceback)




    def run(self,
        seq_fasta_fn,
        complete_traceback):

        # sequences with their ids
        records = parse_fasta(seq_fasta_fn)

        for r in records:
            pass

