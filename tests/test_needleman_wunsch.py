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


def test_guideline_pam():
    """Test cases given on the guideline from 04.02.2019
    """
    nw = NeedlemanWunsch()

    result, info = nw.run("data/xpgma_guideline.fasta",
                    "data/xpgma_guideline.fasta",
                    "data/pam250.txt",
                    False,
                    8,
                    True)

    # the results is a upper triangle matrix of shape n x n.
    seq1_seq2 = result[0][1]
    assert seq1_seq2[3] == 31
    assert len(seq1_seq2[2]) == 1
    assert seq1_seq2[2][0] == ('ILDMDVVEGSAARFDCKVEG_YPDPEVMWFKDDNPVKESRHFQIDYDEEGN', 
                               'RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTD_GRHFVSQTT')

    seq1_seq3 = result[0][2]
    assert seq1_seq3[3] == 44
    assert len(seq1_seq3[2]) == 6
    assert seq1_seq3[2][0] == ('ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN',
                               'ISDTEADIGSNLRWGCAAAGKPRPMVRWLRNGEPL_ASQN_RV__EVLA_')

    seq1_seq4 = result[0][3]
    assert seq1_seq4[3] == 13
    assert len(seq1_seq4[2]) == 24
    assert seq1_seq4[2][0] == ('ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN',
                               'RRLIPAARGGEISILCQPRAAPKATILWSKGTE_ILGNST_RV__TVTSD')

    seq2_seq3 = result[1][2]
    assert seq2_seq3[3] == 15
    assert len(seq2_seq3[2]) == 2
    assert seq2_seq3[2][0] == ('RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT',
                               'ISDTEADIGSNLRWGCAAAGKPRPMV_RWLRNGEP__LASQNR__VEVLA')

    seq2_seq4 = result[1][3]
    assert seq2_seq4[3] == 16
    assert len(seq2_seq4[2]) == 2
    assert seq2_seq4[2][0] == ('RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT',
                               'RRLIPAARGGEISILCQPRAA_PKATILW_SKG_TEILGNSTRVTVT_SD')

    seq3_seq4 = result[2][3]
    assert seq3_seq4[3] == 45
    assert len(seq3_seq4[2]) == 1
    assert seq3_seq4[2][0] == ('ISDTEADIGSNLRWGCAAAGKPRPMVRWLRNGEPLASQNRVEVLA_',
                               'RRLIPAARGGEISILCQPRAAPKATILWSKGTEILGNSTRVTVTSD')


def test_guideline_blosum():
    """Test cases given on the guideline from 04.02.2019
    """
    nw = NeedlemanWunsch()

    result, info = nw.run("data/xpgma_guideline.fasta",
                    "data/xpgma_guideline.fasta",
                    "data/blosum62.txt",
                    False,
                    6,
                    True)

    # the results is a upper triangle matrix of shape n x n.
    seq1_seq2 = result[0][1]
    assert seq1_seq2[3] == 4
    assert len(seq1_seq2[2]) == 8
    assert seq1_seq2[2][0] == ('ILDMDVVEGSAARFDCKVEG_YPDPEVMWFKDDNP__V_KESRHFQIDYDEEGN',
                               'RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHF_V__SQT_T')

    seq1_seq3 = result[0][2]
    assert seq1_seq3[3] == 37
    assert len(seq1_seq3[2]) == 4
    assert seq1_seq3[2][0] == ('ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN',
                               'ISDTEADIGSNLRWGCAAAGKPRPMVRWLRNGEPL_ASQN_RVEV__LA_')

    seq1_seq4 = result[0][3]
    assert seq1_seq4[3] == -4
    assert len(seq1_seq4[2]) == 1
    assert seq1_seq4[2][0] == ('ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN',
                               'RRLIPAARGGEISILCQPRAAPKATILWSKGTEILGNSTRVTVTSD____')

    seq2_seq3 = result[1][2]
    assert seq2_seq3[3] == 3
    assert len(seq2_seq3[2]) == 1
    assert seq2_seq3[2][0] == ('RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT', 
                               'ISDTEADIGSNLRWGC_AAAGKPRPMVRWLRNGEP__LASQNR__VEVLA')


    seq2_seq4 = result[1][3]
    assert seq2_seq4[3] == 9
    assert len(seq2_seq4[2]) == 2
    assert seq2_seq4[2][0] == ('RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT',
                               'RRLIPAARGGEISILCQPRA_APKATILW__SKGTEILGNSTRVTVT_SD')

    seq3_seq4 = result[2][3]
    assert seq3_seq4[3] == 24
    assert len(seq3_seq4[2]) == 1
    assert seq3_seq4[2][0] == ('ISDTEADIGSNLRWGCAAAGKPRPMVRWLRNGEPLASQNRVEVLA_',
                               'RRLIPAARGGEISILCQPRAAPKATILWSKGTEILGNSTRVTVTSD')
    

    