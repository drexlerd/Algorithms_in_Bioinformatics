from prakt.fd import FengDoolittleBase
from xpgma import XPGMA, Node
from needleman_wunsch import NeedlemanWunsch
import math



@FengDoolittleBase.register
class FengDoolittle(FengDoolittleBase):

    def compute_msa(self, 
                    xpgma_root_node : Node,
                    scoring_matrix,
                    cost_gap_open):
        nw = NeedlemanWunsch()
        return self.compute_msa_rec(xpgma_root_node, nw, scoring_matrix, cost_gap_open)


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
                # compute pairwise alignment (similarity)
                score, alignments = nw.compute_optimal_alignments(seq1, seq2, scoring_matrix, cost_gap_open, complete_traceback=False)
                # transform to distance metric

                # select best pairwise alignment according to distance metric

                # make alignment, exchange gap with 'X'

                # return resulting alignment
            else:  # alignment + alignment
                pass


    def similarity_to_distance(self, nw, pairwise_alignment, scoring_matrix, cost_gap_open):
        """Converts similarity score from a pairwise alignment to a distance score
        using approximation algorithm
        
        D(a,b) = - log(S_{a,b}^{eff})

        S_{a,b}^{eff} = (S(a,b) - S_{rand}) / (S_{a,b}^{max} - S_{rand})

        S_{rand} = (1/|A|) * (sum_{x,y in Sigma x Sigma} S(x,y) * N_a(x) * N_b(y)) + gaps(A) * S(-,*)

        S_{a,b}^{max} = (S(a,a) + S(b,b)) / 2
        """
        seq1 = pairwise_alignment[0].replace("_", "")
        seq2 = pairwise_alignment[1].replace("_", "")

        S_ab, _ = nw.compute_optimal_alignments(seq1, seq2, scoring_matrix, cost_gap_open, complete_traceback=False)

        S_aa, _ = nw.compute_optimal_alignments(seq1, seq1, scoring_matrix, cost_gap_open, complete_traceback=False)

        S_bb, _ = nw.compute_optimal_alignments(seq2, seq2, scoring_matrix, cost_gap_open, complete_traceback=False)

        S_ab_max = (S_aa + S_bb) / 2

        S_rand = (1 / len(pairwise_alignment[0])) * \
          sum([scoring_matrix.score(scoring_matrix.alphabet[i], scoring_matrix.alphabet[j]) * self.count_occurences_symbol_in_seq(seq1, scoring_matrix.alphabet[i]) * self.count_occurences_symbol_in_seq(seq2, scoring_matrix.alphabet[j]) for i in range(len(scoring_matrix.alphabet)) for j in range(len(scoring_matrix.alphabet))]) + self.count_gaps_in_pairwise_alignment(pairwise_alignment) * cost_gap_open

        S_eff = (S_ab - S_rand) / (S_ab_max - S_rand)

        return - math.log(S_eff)




    def count_occurences_symbol_in_seq(self, seq, symbol):
        """Counts the number of occurences of symbol in seq

        Args:
          seq (str): A sequence
          symbol (char): A character
        """
        count = 0
        for x in seq:
            if x == symbol:
                count += 1
        return count

    
    def count_gaps_in_pairwise_alignment(self, pairwise_alignment):
        """Counts the number of gaps in the given pairwise alignment. Since Feng-Dootlittle
        exchanges gaps with special character "X" we also need to take care about this special
        character

        Args:
          pairwise_alignment (list(str)): A 2 element list containing a gapped sequences.
                                          Therefore the lengths are equal
                                          
        """
        count = 0
        for i in range(len(pairwise_alignment[0])):
            if pairwise_alignment[0][i] in ["_", "X"] or pairwise_alignment[1][i] in ["_", "X"]:
                count += 1
        return count


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

        xpgma = XPGMA()
        xpgma, n = xpgma.run(seq_fasta_fn, subst_matrix_fn, cost_gap_open, clustering)

