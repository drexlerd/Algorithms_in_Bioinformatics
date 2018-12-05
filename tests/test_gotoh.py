from prakt.gt import GotohBase
from gotoh import Gotoh


def test_instance():
    """Check inheritance."""
    assert issubclass(Gotoh, GotohBase)
    assert isinstance(Gotoh(), GotohBase)


def test_example():
    """Test if run function can be called."""

    gt = Gotoh()
    result = gt.run("data/sequence1.fasta",
                    "data/sequence2.fasta",
                    "data/test_scoring_similarity.txt",
                    False,
                    5,
                    1,  
                    True)

    assert result[0][0][0].id == "idA"
    assert result[0][0][1].id == "idB"
    assert str(result[0][0][0].seq) == "TCCGA"
    assert str(result[0][0][1].seq) == "TACGCAGA"
    assert result[0][0][3] == -3
    assert len(result[0][0][2]) == 3
    assert result[0][0][2][0] == ("T___CCGA",
                             "TACGCAGA")
    assert result[0][0][2][1] == ("TCC___GA",
                             "TACGCAGA")
    assert result[0][0][2][2] == ("TCCG___A",
                             "TACGCAGA")
