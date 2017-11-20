import abc


class NeedlemanWunsch3Base(metaclass=abc.ABCMeta):

    """
    Base class for Needleman-Wunsch implementations.
    """

    @abc.abstractmethod
    def run(self,
            seq_fasta_fn,
            subst_matrix_fn,
            cost_gap_open,
            complete_traceback):
        """
        Calculate optimal alignment(s) with Needleman-Wunsch algorithm for three sequences.

        Args:
            seq_fasta_fn: path to fasta file containing sequences
            subst_matrix_fn: path to substitution matrix
            cost_gap_open: cost to open a gap
            complete_traceback: If True, return all optimal alignments. Otherwise choose a random alignment.

        Returns:
            tuple of
            (id_seq1: fasta id of first sequence,
             seq1: first sequence,
             id_seq2: fasta id of second sequence,
             seq2: second sequence,
             id_seq3: fasta id of second sequence,
             seq3: second sequence,
             score: score of optimal alignment,
             [(aln_string_seq1, aln_string_seq2, aln_string_seq3), ...]: list of tuples containing optimal alignments)
        """
