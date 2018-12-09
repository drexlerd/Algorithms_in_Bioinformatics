from prakt.fd import FengDoolittleBase
from xpgma import XPGMA, Node
from needleman_wunsch import NeedlemanWunsch
from prakt.scoring_func_parser.scoring_func_parser import ScoringMatrix, MetricType
import math



@FengDoolittleBase.register
class FengDoolittle(FengDoolittleBase):

    def compute_msa(self, 
                    xpgma_root_node : Node,
                    nw : NeedlemanWunsch,
                    scoring_matrix,
                    cost_gap_open):
        msa_with_neutral_elements = self.compute_msa_rec(xpgma_root_node, nw, scoring_matrix, cost_gap_open)
        return self.replace_neutral_symbol_with_gap_symbol(msa_with_neutral_elements)
        


    def align_group_to_group_by_group(self, min_i, min_j, group1, group2, group3):
        """Aligns group1 to group2 using group3.

        Args:
          min_i (int): index of the first sequence of group3 in group1
          min_j (int): index of the second sequence of group3 in group2
          group1 (list(str)): An alignment
          group2 (list(str)): An alignment
          group3 (list(str)): An alignment

        Returns:
          list(str): An alignment
        """
        result = ["" for _ in range(len(group1) + len(group2))]
        i1 = 0  # index in alignment 1 best pairwise sequence
        i2 = 0  # index in alignment 2 best pairwise sequence
        for i in range(len(group3[0])):
            c1 = group3[0][i]
            c2 = group3[1][i]
            if c1 == group1[min_i][i1]:
                # copy alignment 1 column to result
                for j in range(len(group1)):
                    result[j] += group1[j][i1]
                i1 += 1
            else:
                for j in range(len(group1)):
                    result[j] += "_"
            if c2 == group2[min_j][i2]:
                for j in range(len(group1), len(group1) + len(group2)):
                    result[j] += group2[j-len(group1)][i2]
                i2 += 1
            else:
                for j in range(len(group1), len(group1) + len(group2)):
                    result[j] += "_"
        return result

        
    def find_best_pairwise_alignment(self, nw, scoring_matrix, cost_gap_open, group1, group2):
        """Compute the best pairwise alignment between sequences of group1, group2.

        Args:
          nw (NeedlemanWunsch): Algorithm for computing pairwise alignments
          group1 (list(str)): An alignment
          group2 (list(str)): An alignment

        Returns:
          list(str): A pairwise alignment
          int: The index of the first sequence of the pairwise alignment in group1
          int: The index of the second sequence of the pairwise alignment in group2
        """
        min_score = float("inf")  # track minimal value
        min_pairwise_alignment = None  # track the best pairwise alignment
        min_i = None
        min_j = None
        for i in range(len(group1)):
            for j in range(len(group2)):
                score, alignments = nw.compute_optimal_alignments(group1[i], group2[j], scoring_matrix, cost_gap_open, complete_traceback=False)
                pairwise_alignment = alignments[0]
                if scoring_matrix.metric_type == MetricType.SIMILARITY:
                    # transform to distance metric                    
                    score = self.similarity_to_distance(nw, pairwise_alignment, scoring_matrix, cost_gap_open)
                if score < min_score:
                    min_score = score
                    min_pairwise_alignment = pairwise_alignment
                    min_i = i
                    min_j = j
        return min_pairwise_alignment, min_i, min_j


    def replace_gap_symbol_with_neutral_symbol(self, group):
        """Replaces all occurences of "_" in the alignment with the neutral symbol "X"

        Args:
          group (list(str)): An alignment
        
        Returns:
          list(str): The alignment with replacement
        """
        return [string.replace("_", "X") for string in group]


    def replace_neutral_symbol_with_gap_symbol(self, group):
        return [string.replace("X", "_") for string in group]


    def compute_msa_rec(self, 
                        node : Node,
                        nw : NeedlemanWunsch,
                        scoring_matrix,
                        cost_gap_open):
        """Recursively computes the multiple sequence alignment.
        Guide tree is traversed depth first.

        Note: If the scoring_matrix is a similarity function, then apply distance transformation

        Base case:
          Leaf node: return alignment consisting of one sequence
        Inductive case:
          Inner node: return alignment between the group of the left child and the group of the right child
                      through the best pairwise alignment

        Args:
          node (Node): The current in traversal
          nw (NeedlemanWunsch): Pairwise alignment algorithm
          scoring_matrix (ScoringMatrix): The scoring matrix object
          cost_gap_open (float): Cost of for a gap 

        Returns:
          list(str): A multiple sequence alignment
        """
        if node.is_leaf():
            group = [str(node.get_seq_record().seq)]
            print("Base case")
            print(group)
            return group  # return the sequence as alignment of 1 sequence
        else:
            assert len(node.get_children()) == 2
            child1, child2 = node.get_children() 
            alignment1 = self.compute_msa_rec(child1.succ, nw, scoring_matrix, cost_gap_open)
            alignment2 = self.compute_msa_rec(child2.succ, nw, scoring_matrix, cost_gap_open)
            print("Inductive case")
            print(alignment1)
            print(alignment2)
            # find pairwise alignment with minimal distance
            min_pairwise_alignment, min_i, min_j = self.find_best_pairwise_alignment(nw, scoring_matrix, cost_gap_open, alignment1, alignment2)
            print(min_pairwise_alignment)
            # align according to best pairwise alignment and return new alignment
            alignment3 = self.align_group_to_group_by_group(min_i, min_j, alignment1, alignment2, min_pairwise_alignment)
            group3 = self.replace_gap_symbol_with_neutral_symbol(alignment3)
            print(alignment3)
            return group3


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
            is_distance_fn,
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
        xpgma, _ = xpgma.run(seq_fasta_fn, subst_matrix_fn, cost_gap_open, clustering)

        nw = NeedlemanWunsch()

        scoring_matrix = ScoringMatrix(subst_matrix_fn, is_distance_fn)

        return self.compute_msa(xpgma, nw, scoring_matrix, cost_gap_open)
