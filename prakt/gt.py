import abc


class GotohBase(metaclass=abc.ABCMeta):

    """
    Base class for Gotoh implementations.
    """

    @abc.abstractmethod
    def run(self,
            seq1_fasta_fn,
            seq2_fasta_fn,
            subst_matrix_fn,
            is_distance_fn,
            affine_cost_gap_open,
            affine_cost_gap_extend,
            complete_traceback):
        """
        Calculate optimal alignment(s) with Needleman-Wunsch algorithm.

        Args:
            seq1_fasta_fn: path to fasta file containing first sequence
            seq2_fasta_fn: path to fasta file containing second sequence
            subst_matrix_fn: path to substitution matrix
            is_distance_fn: True, if a distance matrix is used
            cost_gap_open: cost to open a gap
            cost_gap_extend: cost to extend a gap
            complete_traceback: If True, return all optimal alignments. Otherwise choose a random alignment.

        Returns:
            tuple of
            (id_seq1: fasta id of first sequence,
             seq1: first sequence,
             id_seq2: fasta id of second sequence,
             seq2: second sequence,
             score: score of optimal alignment,
             [(aln_string_seq1, aln_string_seq2), ...]: list of tuples containing optimal alignments)
        """
