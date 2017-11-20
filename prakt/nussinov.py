import abc


class NussinovBase(metaclass=abc.ABCMeta):

    """
    Base class for Nussinov implementations.
    """

    @abc.abstractmethod
    def run(self,
            seq_fasta_fn):
        """
        Fold RNA with Nussinov algorithm.

        Args:
            seq_fasta_fn: path to fasta file containing sequence

        Returns:
            tuple of
            (id_seq: fasta id of sequence,
             seq: sequence,
             structure: dot-bracket string of optimal folding)
        """