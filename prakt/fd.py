import abc


class FengDoolittleBase(metaclass=abc.ABCMeta):

    """
    Base class for Feng-Doolittle implementations.
    """

    @abc.abstractmethod
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