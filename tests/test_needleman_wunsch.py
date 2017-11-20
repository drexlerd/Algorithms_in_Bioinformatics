from prakt.nw import NeedlemanWunschBase
from needleman_wunsch import NeedlemanWunsch


def test_instance():
    """Check inheritance."""
    assert issubclass(NeedlemanWunsch, NeedlemanWunschBase)
    assert isinstance(NeedlemanWunsch(), NeedlemanWunschBase)


def test_example():
    """Example testing the dummy implementation."""

    nw = NeedlemanWunsch()
    result = nw.run("data/sequence1.fa",
                    "data/sequence2.fa",
                    "data/blosum62.txt",
                    5,
                    False)
    (id_seq1, seq1, id_seq2, seq2, score, alignments) = result

    assert id_seq1 == "idA"
    assert id_seq2 == "idB"
    assert seq1 == "FancySequenceA"
    assert seq2 == "FancysequenceB"
    assert score == 1000
    assert alignments[0] == ("Fancy_SequenceA_",
                             "Fancys_equence_B")


# def test_example_fail():
#     """Example of a failing test."""
#
#     nw = NeedlemanWunsch()
#     result = nw.run("data/sequence1.fa",
#                     "data/sequence2.fa",
#                     "data/blosum62.txt",
#                     5,
#                     True)
#     (id_seq1, seq1, id_seq2, seq2, score, alignments) = result
#
#     assert len(alignments) != 1
