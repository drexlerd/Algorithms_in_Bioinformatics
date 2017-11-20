import abc


class XpgmaBase(metaclass=abc.ABCMeta):

    """
    Base class for UPGMA / WPGMA implementations.
    """

    @abc.abstractmethod
    def run(self,
            seq_fasta_fn,
            subst_matrix_fn,
            cost_gap_open,
            clustering):
        """
        Calculate UPGMA / WPGMA clustering.

        Args:
            seq_fasta_fn: path to fasta file containing sequences
            subst_matrix_fn: path to substitution matrix
            cost_gap_open: cost to open a gap
            clustering: select clustering algorithm, either "UPGMA" or "WPGMA"

        Returns:
            NEWICK string with distances
        """