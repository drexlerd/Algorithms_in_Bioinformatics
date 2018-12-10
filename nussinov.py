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
          or c1 == "G" and c2 == "C" or c1 == "C" and c2 == "G" \
          or c1 == "G" and c2 == "U" or c1 == "U" and c2 == "G":
            return True
        return False

    
    def fill_matrix(self, d, sequence, loop_length=1):
        min_loop_length = 2
        for l in range(min_loop_length, len(sequence) + 1):
            for i in range(1, len(sequence) - l + 1):
                j = i + l
                # print("%d, %d" % (i,j))
                # i and j are the indices iterated in diagonal
                # TODO
                case_dist = [(d[i][j-1].value, (i, j, None), Case.UNPAIRED)]
                
                for k in range(i, j - loop_length):
                    if self.is_base_pair(sequence[k-1],sequence[j-1]):
                        # print("%d, %d, %s, %s" % (i, j, sequence[k-1], sequence[j-1]))
                        case_dist.append((d[i][k-1].value + d[k+1][j-1].value + 1, (i, j, k), Case.PAIRED_SUCCESS))
                    else:
                        case_dist.append((0, (i, j, k), Case.PAIRED_FAIL))

                result_sequence = Max(case_dist)

                # store result of the recursion
                d[i][j].SetValue(result_sequence[0][0])  # set value
                # set predecessors
                for result_case in result_sequence:
                    d[i][j].AddPredecessor(result_case[1], result_case[2])
        
        #print()
        #for i in range(len(sequence)):
        #    for j in range(len(sequence) + 1):
        #        print("%5.2f" % (d[i][j].value), end='')
        #    print()


    def traceback_rec(self, d, i, j, structure_index=0, structures=[(0, None)]):
        """Initial call: traceback_rec(d, 0, len(sequence), 0, [(0, None)])

        Args:
          d (list(list(Cell))): DP matrix (with traceback pointers)
          i (int): smallest index in sequence (included)
          j (int): biggest index in sequence (included)
          structure_index (int): the index in structures we are currently at in our recursion
          structure (list(tuple(int,tuple(int,int)))): pointer to head recursion + a tuples of base pairings
        """
        current_cell = d[i][j]
        if not self._is_traceback_base_case(current_cell):
            for (i, j, k), case in current_cell.pre:
                # copy the structure because of multiple predecessor
                if case == Case.UNPAIRED:
                    self.traceback_rec(d, i, j-1, structure_index, structures)
                elif case == Case.PAIRED_SUCCESS:
                    structures.append((structure_index, (j, k)))
                    self.traceback_rec(d, i, k-1, len(structures) - 1, structures)
                    self.traceback_rec(d, k + 1, j - 1, len(structures) - 1, structures)
                elif case == Case.PAIRED_FAIL:
                    self.traceback_rec(d, i, k-1, len(structures) - 1, structures)
                    self.traceback_rec(d, k + 1, j - 1, len(structures) - 1, structures)
            return structures


    def _is_traceback_base_case(self, cell):
        if cell.i + 2 > cell.j:
            return True
        return False


    def compute_optimal_abstract_structure(self, sequence, complete_traceback):
        d = [[Cell(i, j) for j in range(len(sequence) + 1)] for i in range(len(sequence) + 1)]

        self.fill_matrix(d, sequence, loop_length=1)

        abstract_structures = self.traceback_rec(d, 1, len(sequence), 0, [(0, None)]) 

        return abstract_structures, d[1][len(sequence)].value
        

    def convert_abstract_structure_to_structure(self, sequence, amount_pairs, abstract_structure):
        structures = []
        # track which abstract structures are already consumed
        for index, (j, k) in reversed(abstract_structure[1:]):
            count = 0
            r = ["."] * len(sequence)
            while True:
                r[j-1] = ")"
                r[k-1] = "("
                count += 1
                if count == amount_pairs:
                    structures.append("".join(r))
                if abstract_structure[index][1] == None:
                    break
                index, (j, k) = abstract_structure[index]
        return structures





    def run(self,
        seq_fasta_fn,
        complete_traceback):

        # sequences with their ids
        records = parse_fasta(seq_fasta_fn)

        for r in records:
            pass

