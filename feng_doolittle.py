from prakt.fd import FengDoolittleBase
from xpgma import XPGMA, Node
from needleman_wunsch import NeedlemanWunsch




@FengDoolittleBase.register
class Gotoh(FengDoolittleBase):

    def compute_msa(self, 
                    xpgma_root_node : Node,
                    scoring_matrix,
                    cost_gap_open):
        nw = NeedlemanWunsch()
        return self.compute_msa(xpgma_root_node, nw, scoring_matrix, cost_gap_open)


    def compute_msa_rec(self, 
                        node : Node,
                        nw : NeedlemanWunsch,
                        scoring_matrix,
                        cost_gap_open):
        if node.is_leaf():
            return node.get_seq_record()
        else:
            assert len(node.get_children) == 2
            child1, child2 = node.get_children()
            if child1.succ.is_leaf() or child2.succ.is_leaf():  # alignment + sequence
                pass
            elif child1.succ.is_leaf and child2.succ.is_leaf():  # sequence + sequence
                seq1_record = child1.get_seq_record()
                seq2_record = child2.get_seq_record()
                seq1 = str(seq1_record.seq)
                seq2 = str(seq2_record.seq)

                score, alignments = nw.compute_optimal_alignments(seq1, seq2, scoring_matrix, cost_gap_open, complete_traceback=False)
            else:  # alignment + alignment
                pass


    def run(self,
            seq_fasta_fn,
            subst_matrix_fn,
            cost_gap_open,
            clustering):
        """
        Calculate optimal alignment with Feng-Doolittle algorithm.

        Args:
            seq_fasta_fn: path to fasta file containing sequences
            subst_matrix_fn: path to substitution matrix
            cost_gap_open: cost to open a gap
            clustering: select clustering algorithm, either "UPGMA" or "WPGMA"

        Returns:
            tuple of
            (score: sum-of-pairs score of optimal alignment,
             [aln_seq1, aln_seq2, ...]: final alignment as list of strings
             [aln_seq1_id, aln_seq2_id, ...]: list of sequence ids in same order as aligned sequences)
        """

        xpgma = XPGMA
        xpgma, n = xpgma.run(seq_fasta_fn, subst_matrix_fn, cost_gap_open, clustering)