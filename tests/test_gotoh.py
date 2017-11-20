from prakt.gt import GotohBase
from gotoh import Gotoh


def test_instance():
    """Check inheritance."""
    assert issubclass(Gotoh, GotohBase)
    assert isinstance(Gotoh(), GotohBase)


def test_example():
    """Test if run function can be called."""

    gt = Gotoh()
    result = gt.run("data/sequence1.fa",
                    "data/sequence2.fa",
                    "data/blosum62.txt",
                    5,
                    20,
                    False)
    (id_seq1, seq1, id_seq2, seq2, score, alignments) = result

    assert id_seq1 == "idA"
    assert id_seq2 == "idB"
    assert seq1 == "FancySequenceA"
    assert seq2 == "FancysequenceB"
    assert score == 1000
    assert alignments[0] == ("Fancy_SequenceA_",
                             "Fancys_equence_B")
