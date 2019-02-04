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
    assert result[0][0][3] == -4
    assert len(result[0][0][2]) == 3
    assert result[0][0][2][0] == ('T___CCGA', 'TACGCAGA') 


def test_guideline_pam():
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
    assert seq1_seq2[3] == 33
    assert len(seq1_seq2[2]) == 2
    assert seq1_seq2[2][0] == ('ILDMDVVEGSAARFDCKVEG_YPDPEVMWFKDDNP___VKESRHFQIDYDEEGN',
                               'RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHF____VSQTT')

    seq1_seq3 = result[0][2]
    assert seq1_seq3[3] == 60
    assert len(seq1_seq3[2]) == 2
    assert seq1_seq3[2][0] == ('ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN',
                               'ISDTEADIGSNLRWGCAAAGKPRPMVRWLRNGEPLASQNRVEV_____LA')

    seq1_seq4 = result[0][3]
    assert seq1_seq4[3] == 30
    assert len(seq1_seq4[2]) == 2
    assert seq1_seq4[2][0] == ('ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN',
                               'RRLIPAARGGEISILCQPRAAPKATILWSKGTEILGNSTRVTV____TSD')

    seq2_seq3 = result[1][2]
    assert seq2_seq3[3] == 17
    assert len(seq2_seq3[2]) == 2
    assert seq2_seq3[2][0] == ('RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT',
                               'ISDTEADIGSNLRWGCAAAGKPRPMV_RWLRNG____EPLASQNRVEVLA')

    seq2_seq4 = result[1][3]
    assert seq2_seq4[3] == 9
    assert len(seq2_seq4[2]) == 1
    assert seq2_seq4[2][0] == ('RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT',
                               'RRLIPAARGGEISILCQPRAA_PKATILW__SKGTEILGNSTRVTVT_SD')

    seq3_seq4 = result[2][3]
    assert seq3_seq4[3] == 41
    assert len(seq3_seq4[2]) == 1
    assert seq3_seq4[2][0] == ('ISDTEADIGSNLRWGCAAAGKPRPMVRWLRNGEPLASQNRVEVLA_',
                               'RRLIPAARGGEISILCQPRAAPKATILWSKGTEILGNSTRVTVTSD')


def test_guideline_blosum():
    """Test cases given on the guideline from 04.02.2019
    
    the results are not as provided by the guideline
    """
    gotoh = Gotoh()

    result, info = gotoh.run("data/xpgma_guideline.fasta",
                    "data/xpgma_guideline.fasta",
                    "data/blosum62.txt",
                    False,
                    11,
                    1,
                    True)

    # the results is a upper triangle matrix of shape n x n.
    seq1_seq2 = result[0][1]
    assert seq1_seq2[3] == 0
    assert len(seq1_seq2[2]) == 1
    assert seq1_seq2[2][0] == ('ILDMDVVEGSAARFDCKVEG_YPDPEVMWFKDDNP___VKESRHFQIDYDEEGN',
                               'RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFV____SQTT')

    seq1_seq3 = result[0][2]
    assert seq1_seq3[3] == 41
    assert len(seq1_seq3[2]) == 3
    assert seq1_seq3[2][0] == ('ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN',
                               'ISDTEADIGSNLRWGCAAAGKPRPMVRWLRNGEPLASQNRVEV_____LA')

    seq1_seq4 = result[0][3]
    assert seq1_seq4[3] == 5
    assert len(seq1_seq4[2]) == 1
    assert seq1_seq4[2][0] == ('ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN',
                               'RRLIPAARGGEISILCQPRAAPKATILWSKGTEILGNSTRVTVTSD____')

    seq2_seq3 = result[1][2]
    assert seq2_seq3[3] == -5
    assert len(seq2_seq3[2]) == 4
    assert seq2_seq3[2][0] == ('RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT',
                               'ISDTEADIGSNLRWGC_AAAGKPRPMVRWLRNG____EPLASQNRVEVLA')

    seq2_seq4 = result[1][3]
    assert seq2_seq4[3] == -4
    assert len(seq2_seq4[2]) == 2
    assert seq2_seq4[2][0] == ('RDPVKTHEGWGVMLPCNPPAHYPGLSYRWLLNEFPNFIPTDGRHFVSQTT',
                               'RRLIPAARGGEISILCQPRA_APKATILW__SKGTEILGNSTRVTVT_SD')

    seq3_seq4 = result[2][3]
    assert seq3_seq4[3] == 18
    assert len(seq3_seq4[2]) == 1
    assert seq3_seq4[2][0] == ('ISDTEADIGSNLRWGCAAAGKPRPMVRWLRNGEPLASQNRVEVLA_',
                               'RRLIPAARGGEISILCQPRAAPKATILWSKGTEILGNSTRVTVTSD')