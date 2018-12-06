from prakt.nw import XpgmaBase
from needleman_wunsch import NeedlemanWunsch
from prakt.fasta_parser.fasta_parser import parse_fasta
from prakt.scoring_func_parser.scoring_func_parser import ScoringMatrix


class Cell(object):
    def __init__(self, i, j, distance):
        self.distance = 0
        self.succ = None 
        self.i = i
        self.j = j


class Node(object):
    def __init__(self, seq_record):
        self.seq_record = seq_record
        self.children = None


class Edge(object):
    def __init__(self):
        self.weight = 0
        self.succ = None
            

@XpgmaBase.register
class XPGMA(XpgmaBase):

    def generate_upgma(self, m, l, n):
        """
        Args:
          m (list(list(int))): cluster distance matrix
          n (dict(int, Node)): Cluster distance matrix index to node mapping
          l (list(int)): iterationlist
        """
        while len(l) > 1:  # stop if only 1 cluster remains
            pass

    def run(self,
            seq_fasta_fn,
            subst_matrix_fn,
            cost_gap_open,
            clustering):
            scoring_matrix = ScoringMatrix(subst_matrix_fn, True)

            seq_records = parse_fasta(seq_fasta_fn)
            seqs = [str(x.seq) for x in seq_records]

            # cluster distance matrix, containing pairwise distance information
            m_size = 2 * len(seqs) - 1  # additional len(seqs) - 1 rows when merging clusters
            m = [[0 for i in range(m_size)] for j in range(m_size)]

            # iterationlist, containing the matrix row/col indices of the current clusters
            # this is used to avoid having to clean the matrix after merge
            l = [i for i in range(len(seqs))]  # initially only singleton clusters

            # cluster distance matrix index to Node mapping
            initial_cluster = [Node(seq_records[i]) for i in range(len(seq_records))]
            n = dict(zip(list(range(len(initial_cluster)), initial_cluster)))

            # compute pairwise distances using NW
            # Note: no check if matrix is distance matrix
            nw = NeedlemanWunsch()

            result, info = nw.run(seq_fasta_fn,
                    seq_fasta_fn,
                    subst_matrix_fn,
                    True,
                    cost_gap_open,
                    False)

            # initialize cluster distance matrix with computed distances
            for i in range(len(seqs)):
                for j in range(len(seqs)):
                    m[i][j] = result[i][j][3]


            

            




            
