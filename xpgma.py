from prakt.xpgma import XpgmaBase
from needleman_wunsch import NeedlemanWunsch
from prakt.fasta_parser.fasta_parser import parse_fasta
from prakt.scoring_func_parser.scoring_func_parser import ScoringMatrix


class Node(object):
    def __init__(self, seq_record=None):
        self.seq_record = seq_record
        self.children = None
        self.distance = 0  # distance to leaf node (Note: this is enough for UPGMA but not for WPGMA)
    
    def set_children(self, children_list):
        self.children = children_list

    def set_distance(self, distance):
        self.distance = distance
        
    def get_distance(self):
        return self.distance

    def __repr__(self):
        return "distance=%.2f\n" % (self.distance)


class Edge(object):
    def __init__(self, weight=0, succ=None):
        self.weight = weight
        self.succ = succ
            

@XpgmaBase.register
class XPGMA(XpgmaBase):
    def find_smallest_distance_clusters(self, m, l):
        """Finds the two clusters which have minimal distance
        """
        min_distance = float("inf")
        min_ci = None
        min_cj = None
        print(l)
        for i in range(len(l)):
            for j in range(i+1, len(l)):
                c1 = l[i]
                c2 = l[j]
                if m[c1][c2] <= min_distance:
                    min_distance = m[c1][c2]
                    min_ci = c1
                    min_cj = c2
        return min_ci, min_cj, min_distance

                
    def generate_upgma(self, m, l, n):
        """
        Args:
          m (list(list(int))): cluster distance matrix (upper triangle matrix)
          n (dict(int, Node)): Cluster distance matrix index to node mapping
          l (list(int)): iterationlist
        """
        for new_cluster_index in range(len(l), len(m)):  # compute all the cluster merges
            ci, cj, min_distance = self.find_smallest_distance_clusters(m, l)
            for i in range(len(l)):
                for j in range(len(l)):
                    print("%5.1f" % (m[l[i]][l[j]]), end="")
                print()
            # now merge ci and cj
            new_cluster_node = Node()
            # set edge weights (according to UPGMA)
            e1 = Edge(weight=min_distance/2 - n[ci].get_distance(), succ=n[ci])
            e2 = Edge(weight=min_distance/2 - n[cj].get_distance(), succ=n[cj])
            # set child edges
            new_cluster_node.set_children([e1, e2])
            new_cluster_node.set_distance(min_distance/2)
            # add node in dict
            n[new_cluster_index] = new_cluster_node
            # remove indices from l
            l.remove(ci)
            l.remove(cj)
            # set distances to other clusters (according to UPGMA)
            for ck in l:
                # need to take max to end up with the value in upper triangle matrix
                m[ck][new_cluster_index] = (max(m[ck][ci], m[ci][ck]) + max(m[ck][cj], m[cj][ck])) / 2
            # add new index
            l.append(new_cluster_index)
        return new_cluster_node, n



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
                for j in range(i+1, len(seqs)):
                    m[i][j] = result[i][j][3]


            

            




            
