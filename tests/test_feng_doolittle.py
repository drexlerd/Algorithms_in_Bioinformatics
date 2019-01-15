from prakt.fd import FengDoolittleBase
from feng_doolittle import FengDoolittle
from needleman_wunsch import NeedlemanWunsch
from prakt.scoring_func_parser.scoring_func_parser import ScoringMatrix
from prakt.util.util import similarity_to_distance, count_occurences_symbol_in_seq, count_gaps_in_pairwise_alignment
import math


def test_count_occurences_symbol_in_seq():
    seq1 = ("__TCCGA_", "TACGCAGA")[0]
    seq2 = ("__TCCGA_", "TACGCAGA")[1]
    fd = FengDoolittle()
    count = count_occurences_symbol_in_seq(seq1, "C")
    assert count == 2
    count = count_occurences_symbol_in_seq(seq1, "G")
    assert count == 1


def test_similarity_to_distance():
    scoring_matrix = ScoringMatrix("data/test_scoring_similarity.txt", is_distance_fn=False, cost_gap_open=1)
    nw = NeedlemanWunsch()
    fd = FengDoolittle()
    pairwise_alignment = ("__TCCGA_", "TACGCAGA")
    distance = similarity_to_distance(nw, pairwise_alignment, scoring_matrix)

    count = count_gaps_in_pairwise_alignment(pairwise_alignment)
    assert count == 3

    # the right hand side was computed from hand and is - log(S_eff)
    assert distance == - math.log((2 -- 14/8) / (6.5 -- 14/8))

def test_feng_doolittle():
    fd = FengDoolittle()

    # test with low gap cost
    msa = fd.run("data/xpgma.fasta",
            "data/test_scoring_distance.txt",
            True,
            1,
            "wpgma")
    assert msa == ['GCT____TGTTACGAT', 'TC_____TGTTACGAT', 'ACTTGACCG_TT___T', 'ACTACACCCTTATGAG', 'ACTTGTCCGAAACGAT', 'AGATGACCGTTTCGAT']
    
    # test with high gap cost => likelier to align gap with gap
    msa = fd.run("data/xpgma.fasta",
            "data/test_scoring_distance.txt",
            True,
            2,
            "wpgma")
    
    assert msa == ['ACTACACCCTTATGAG', 'ACTTGTCCGAAACGAT', 'AGATGACCGTTTCGAT', 'ACT____TGACCGTTT', 'GCT____TGTTACGAT', 'TC_____TGTTACGAT']

    
def test_feng_doolittle_similarity():
    fd = FengDoolittle()

    msa = fd.run("data/xpgma.fasta",
            "data/test_scoring_similarity.txt",
            False,
            1,
            "wpgma")

def test_feng_doolittle_guideline():
	fd = FengDoolittle()

	msa = fd.run("data/feng_doolittle_guideline.fasta",
				 "data/blosum62.txt",
				 False,
				 8,
				 "wpgma")

	