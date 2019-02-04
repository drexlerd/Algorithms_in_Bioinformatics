from prakt.nw import NeedlemanWunschBase
from needleman_wunsch import NeedlemanWunsch


def test_instance():
    """Check inheritance."""
    assert issubclass(NeedlemanWunsch, NeedlemanWunschBase)
    assert isinstance(NeedlemanWunsch(), NeedlemanWunschBase)


def test_example_distance():
    """Test using distance scoring function"""

    nw = NeedlemanWunsch()
    result, info = nw.run("data/sequence1.fasta",
                    "data/sequence2.fasta",
                    "data/test_scoring_distance.txt",
                    True,
                    1,
                    True)

    assert result[0][0][0].id == "idA"
    assert result[0][0][1].id == "idB"
    assert str(result[0][0][0].seq) == "TCCGA"
    assert str(result[0][0][1].seq) == "TACGCAGA"
    assert result[0][0][3] == -2
    assert len(result[0][0][2]) == 1
    assert result[0][0][2][0] == ("T_C_C_GA",
                             "TACGCAGA")


def test_example_similarity():
    """Test using similarity scoring function
    """

    nw = NeedlemanWunsch()

    result, info = nw.run("data/sequence1.fasta",
                    "data/sequence2.fasta",
                    "data/test_scoring_similarity.txt",
                    True,
                    1,
                    True)

    assert result[0][0][0].id == "idA"
    assert result[0][0][1].id == "idB"
    assert str(result[0][0][0].seq) == "TCCGA"
    assert str(result[0][0][1].seq) == "TACGCAGA"
    assert result[0][0][3] == 4
    assert len(result[0][0][2]) == 6
    assert result[0][0][2][0] == ("__TCCGA_",
                                  "TACGCAGA")

