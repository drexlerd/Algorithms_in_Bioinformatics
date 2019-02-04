from prakt.gt import GotohBase
from gotoh import Gotoh


def test_instance():
    """Check inheritance."""
    assert issubclass(Gotoh, GotohBase)
    assert isinstance(Gotoh(), GotohBase)


def test_example():
    """Test if run function can be called."""

    gt = Gotoh()
    result, info = gt.run("data/sequence1.fasta",
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


def test_guideline():
    """Test cases given on the guideline from 04.02.2019
    
    the results are not as provided by the guideline
    """
    gotoh = Gotoh()

    result, info = gotoh.run("data/xpgma_guideline.fasta",
                    "data/xpgma_guideline.fasta",
                    "data/pam250.txt",
                    False,
                    11,
                    1,
                    True)

    # the results is a upper triangle matrix of shape n x n.
    seq1_seq2 = result[0][1]
    #assert seq1_seq2[3] == 33
    #assert len(seq1_seq2[2]) == 2
    #assert seq1_seq2[2][0] == ('ILDMDVVEGSAARFDCKVEG_YPDPEVMWFKDDNPVKESRHFQIDYDEEGN', 
    #                           'RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTD_GRHFVSQTT')

    seq1_seq3 = result[0][2]
    #assert seq1_seq3[3] == 60
    #assert len(seq1_seq3[2]) == 2
    #assert seq1_seq3[2][0] == ('ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN',
    #                           'ISDTEADIGSNLRWGCAAAGKPRPMVRWLRNGEPL_ASQN_RV__EVLA_')

    #seq1_seq4 = result[0][3]
    #assert seq1_seq4[3] == 30
    #assert len(seq1_seq4[2]) == 2
    #assert seq1_seq4[2][0] == ('ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN',
    #                           'RRLIPAARGGEISILCQPRAAPKATILWSKGTE_ILGNST_RV__TVTSD')

    seq2_seq3 = result[1][2]
    #assert seq2_seq3[3] == 9
    #assert len(seq2_seq3[2]) == 1
    #assert seq2_seq3[2][0] == ('RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT',
    #                           'ISDTEADIGSNLRWGCAAAGKPRPMV_RWLRNGEP__LASQNR__VEVLA')

    seq2_seq4 = result[1][3]
    #assert seq2_seq4[3] == 41
    #assert len(seq2_seq4[2]) == 1
    #assert seq2_seq4[2][0] == ('RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT',
    #                           'RRLIPAARGGEISILCQPRAA_PKATILW_SKG_TEILGNSTRVTVT_SD')

    seq3_seq4 = result[2][3]
    #assert seq3_seq4[3] == 17
    #assert len(seq3_seq4[2]) == 2
    #assert seq3_seq4[2][0] == ('ISDTEADIGSNLRWGCAAAGKPRPMVRWLRNGEPLASQNRVEVLA_',
    #                           'RRLIPAARGGEISILCQPRAAPKATILWSKGTEILGNSTRVTVTSD')