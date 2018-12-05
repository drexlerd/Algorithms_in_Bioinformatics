"""
author: dominik drexler <drexlerd@informatik.uni-freiburg.de>
"""

from Bio import SeqIO
from enum import Enum

PROTEIN_ALPHABET = ['A', 'R', 'N', 'D', 'B', 'C', 'E', 'Q', 'Z', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
DNA_ALPHABET = ['A', 'G', 'C', 'T']

class SequenceType(Enum):
    DNA = 0,
    PROTEIN = 1,


def check_sequences_alphabet(sequence_list, sequence_type):
    """
    Args:
      sequence_list (list(SeqIO)): A list of SeqIO objects
      sequence_type (SequenceType): The type of the sequence to check against
    """
    for record in sequence_list:
        for c in str(record.seq):
            if sequence_type == SequenceType.PROTEIN:
                if c not in PROTEIN_ALPHABET:
                    return False
            elif sequence_type == SequenceType.DNA:
                if c not in DNA_ALPHABET:
                    return False
    return True



def parse_fasta(filename):
    """For now return just the first sequence
    """

    return list(SeqIO.parse(filename, "fasta")) 
