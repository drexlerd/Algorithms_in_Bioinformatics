from prakt.nw import XpgmaBase
from needleman_wunsch import NeedlemanWunsch
from prakt.fasta_parser.fasta_parser import parse_fasta
from prakt.scoring_func_parser.scoring_func_parser import ScoringMatrix


@XpgmaBase.register
class XPGMA(XpgmaBase):
    def find_cell_indices_with_lowest_pairwise_distance(self, matrix):
        return min([matrix[i][j] for i in range(len(matrix)) for j in range(len(matrix))], key=x : x.distance)


    def merge_sequences_from_cell(self, cell, matrix):
        pass


    def run(self,
            seq_fasta_fn,
            subst_matrix_fn,
            cost_gap_open,
            clustering):
            scoring_matrix = ScoringMatrix(subst_matrix_fn, True)

            seq_records = parse_fasta(seq1_fasta_fn)
            seqs = [str(x.seq) for x in seq_records]

            M = [[0 for i in range(len(seqs))] for j in range(len(seqs))]

            nw = NeedlemanWunsch()  # for computing distances

            result = nw.run(seq_fasta_fn,
                    seq_fasta_fn,
                    subst_matrix_fn,
                    True,
                    cost_gap_open,
                    False)

            self.m = [[Cell(i, j, result[i][j][3]) for j in range(len(seqs))] for i in range(len(seqs))]

            

class Cell(object):
    def __init__(self, i, j, distance):
        self.distance = 0
        self.succ = None 
        self.i = i
        self.j = j


class Node(object):
    def __init__(self):
        self.children = None


class Edge(object):
    def __init__(self):
        self.weight = 0
        self.succ = None
            
