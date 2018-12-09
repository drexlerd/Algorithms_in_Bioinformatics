from prakt.fd import FengDoolittleBase
from feng_doolittle import FengDoolittle
from needleman_wunsch import NeedlemanWunsch
from prakt.scoring_func_parser.scoring_func_parser import ScoringMatrix
import math


def test_count_occurences_symbol_in_seq():
    seq1 = ("__TCCGA_", "TACGCAGA")[0]
    seq2 = ("__TCCGA_", "TACGCAGA")[1]
    fd = FengDoolittle()
    count = fd.count_occurences_symbol_in_seq(seq1, "C")
    assert count == 2
    count = fd.count_occurences_symbol_in_seq(seq1, "G")
    assert count == 1


def test_similarity_to_distance():
    scoring_matrix = ScoringMatrix("data/test_scoring_similarity.txt", is_distance_fn=False)
    nw = NeedlemanWunsch()
    fd = FengDoolittle()
    pairwise_alignment = ("__TCCGA_", "TACGCAGA")
    cost_gap_open = -1  # in the similarity mesasure this gets converted from 1 to -1
    distance = fd.similarity_to_distance(nw, pairwise_alignment, scoring_matrix, cost_gap_open)

    count = fd.count_gaps_in_pairwise_alignment(pairwise_alignment)
    assert count == 3

    # the right hand side was computed from hand and is - log(S_eff)
    assert distance == - math.log((2 -- 14/8) / (6.5 -- 14/8))

def test_feng_doolittle():
    fd = FengDoolittle()

    msa = fd.run("data/xpgma.fasta",
            "data/test_scoring_distance.txt",
            True,
            1,
            "wpgma")

    print(msa)