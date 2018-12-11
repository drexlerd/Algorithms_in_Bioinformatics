"""
Module name: XPGMA
Module author: dominik drexler <drexlerd@informatik.uni-freiburg.de>
"""

from prakt.xpgma import XpgmaBase
from needleman_wunsch import NeedlemanWunsch
from prakt.fasta_parser.fasta_parser import parse_fasta
from prakt.scoring_func_parser.scoring_func_parser import ScoringMatrix, MetricType
from prakt.util.util import similarity_to_distance
import argparse


class Node(object):
    """This class represents a leaf node or an inner node in a XPGMA
    """
    def __init__(self, seq_record=None):
        # leaf nodes contain a sequence record
        self.seq_record = seq_record
        # children for traversing
        self.children = None
        # distance for setting weights of new edges 
        self.distance = 0
        # cluster size needed in upgma computation
        self.cluster_size = 1
    
    def set_children(self, children_list):
        self.children = children_list

    def set_distance(self, distance):
        self.distance = distance

    def set_cluster_size(self, value):
        self.cluster_size = value

    def get_cluster_size(self):
        return self.cluster_size
        
    def get_distance(self):
        return self.distance
        
    def get_seq_record(self):
        return self.seq_record

    def get_children(self):
        return self.children
        
    def is_leaf(self):
        if self.children is None:
            return True
        return False
        
    def __repr__(self):
        return "distance=%.2f\n" % (self.distance)


class Edge(object):
    """This class represents an edge in a XPGMA
    """
    def __init__(self, weight=0, succ=None):
        self.weight = weight
        self.succ = succ

    def __repr__(self):
        return "Edge weight: %.2f" % (self.weight)
            

@XpgmaBase.register
class XPGMA(XpgmaBase):
    def find_smallest_distance_clusters(self, m, l):
        """Finds the two clusters which have minimal distance

        Args:
          m (list(list(int))): cluster distance matrix (upper triangle matrix)
          l (list(int)): list containing the cluster indices of the root clusters
        """
        min_distance = float("inf")
        min_ci = None
        min_cj = None
        # print(l)
        for i in range(len(l)):
            for j in range(i+1, len(l)):
                c1 = l[i]
                c2 = l[j]
                if m[c1][c2] <= min_distance:
                    min_distance = m[c1][c2]
                    min_ci = c1
                    min_cj = c2
        return min_ci, min_cj, min_distance

                
    def generate_wpgma(self, m, l, n):
        """Compute a WPGMA

        Args:
          m (list(list(int))): cluster distance matrix (upper triangle matrix)
          n (dict(int, Node)): cluster index to node mapping
          l (list(int)): list containing the cluster indices of the root clusters

        Returns:
          new_cluster_node (Node): Root node of the XPGMA
          n
        """
        for new_cluster_index in range(len(l), len(m)):  # compute all the cluster merges
            ci, cj, min_distance = self.find_smallest_distance_clusters(m, l)
            #for i in range(len(l)):
            #    for j in range(len(l)):
            #        print("%5.1f" % (m[l[i]][l[j]]), end="")
            #    print()
            # now merge ci and cj
            new_cluster_node = Node()
            # new_cluster_node.set_cluster_size(n[ci].get_cluster_size() + n[cj].get_cluster_size())
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
                m[ck][new_cluster_index] = ((m[ck][ci] + m[ci][ck]) + (m[ck][cj] + m[cj][ck])) / 2
            # add new index
            l.append(new_cluster_index)
        return new_cluster_node, n


    def generate_upgma(self, m, l, n):
        """Computes a UPGMA

        Args:
          m (list(list(int))): cluster distance matrix (upper triangle matrix)
          n (dict(int, Node)): cluster index to node mapping
          l (list(int)): list containing the cluster indices of the root clusters

        Returns:
          new_cluster_node (Node): Root node of the XPGMA
          n
        """
        for new_cluster_index in range(len(l), len(m)):  # compute all the cluster merges
            ci, cj, min_distance = self.find_smallest_distance_clusters(m, l)
            #for i in range(len(l)):
            #    for j in range(len(l)):
            #        print("%5.1f" % (m[l[i]][l[j]]), end="")
            #    print()
            # merge ci and cj
            new_cluster_node = Node()
            ci_cluster_size = n[ci].get_cluster_size()
            cj_cluster_size = n[cj].get_cluster_size()
            new_cluster_node.set_cluster_size(ci_cluster_size + cj_cluster_size)
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
                m[ck][new_cluster_index] = (ci_cluster_size * (m[ck][ci] + m[ci][ck]) + cj_cluster_size * (m[ck][cj] + m[cj][ck])) / (ci_cluster_size + cj_cluster_size)
            # add new index
            l.append(new_cluster_index)
        return new_cluster_node, n


    def run(self,
            seq_fasta_fn,
            subst_matrix_fn,
            is_distance_fn,
            cost_gap_open,
            clustering):
            """
            Computes a XPGMA

            Args:
              seq_fasta_fn (str): The relative path to a fasta file
              subst_matrix_fn (str): The relative path to a scoring matrix file
              is_distance_fn (bool): If True, handle scoring matrix as distance measure, else similarity measure
              cost_gap_open (int): gap cost open
              clustering (str): either "upgma" or "wpgma"

            Returns:
                new_cluster_node (Node): Root node of the XPGMA
                n
            """
            scoring_matrix = ScoringMatrix(subst_matrix_fn, is_distance_fn, cost_gap_open)
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
            n = dict(zip(list(range(len(initial_cluster))), initial_cluster))

            # compute pairwise distances using NW
            # Note: no check if matrix is distance matrix
            nw = NeedlemanWunsch()

            result, info = nw.run(seq_fasta_fn,
                    seq_fasta_fn,
                    subst_matrix_fn,
                    is_distance_fn,
                    cost_gap_open,
                    False)

            # initialize cluster distance matrix with computed distances
            for i in range(len(seqs)):
                for j in range(i+1, len(seqs)):
                    if scoring_matrix.metric_type == MetricType.DISTANCE:
                        m[i][j] = result[i][j][3]
                    else:
                        m[i][j] == similarity_to_distance(nw, result[i][j][2][0], scoring_matrix)

            if clustering == "wpgma":
                return self.generate_wpgma(m, l, n)
            elif clustering == "upgma":
                return self.generate_upgma(m, l, n)


            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="XPGMA command line tool")
    parser.add_argument("seq_fasta_fn", type=str)
    parser.add_argument("subst_matrix_fn", type=str)
    parser.add_argument("cost_gap_open", type=int)
    parser.add_argument('clustering', choices=['wpgma', 'upgma'])
    parser.add_argument("--d", "--is_distance_fn", action='store_true')
    args = parser.parse_args()

    xpgma = XPGMA()

    xpgma, n = xpgma.run(args.seq_fasta_fn, args.subst_matrix_fn, args.is_distance_fn, args.cost_gap_open, args.clustering)

    print(xpgma)

            




            
