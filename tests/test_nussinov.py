from nussinov import Nussinov
from prakt.basis_classes.cell import Cell

def test_fill_matrix():
    nussinov = Nussinov()

    sequence = "GGUCCAC"

    d = [[Cell(i, j) for j in range(len(sequence) + 1)] for i in range(len(sequence) + 1)]

    nussinov.fill_matrix(d, sequence)


def test_nussinov_guideline():
    nussinov = Nussinov()

    structures, amount_pairs = nussinov.run("data/nussinov_guideline.fasta", True)
    
    structures_seq1 = structures[0]
    amount_pair_seq1 = amount_pairs[0]
    assert amount_pair_seq1 == 11
    assert '(((((()())(.((()..))..)...))))' in structures_seq1

    structures_seq2 = structures[1]
    amount_pair_seq2 = amount_pairs[1]
    assert amount_pair_seq2 == 24
    assert '((((.((((((((((...)))))..)(((()))))()()))((((...)))))...))))' in structures_seq2

    structures_seq3 = structures[2]
    amount_pair_seq3 = amount_pairs[2]
    assert amount_pair_seq3 == 26
    assert '(..(((((.())())((.((((()))(((()())))(().))))(.((())))).).))).' in structures_seq3

    structures_seq4 = structures[3]
    amount_pair_seq4 = amount_pairs[3]
    assert amount_pair_seq4 == 27
    assert '(((((()(()(((.....(.(((((.()((.)))((()()))).)((().)).)..))))))).))))' in structures_seq4

    structures_seq5 = structures[4]
    amount_pair_seq5 = amount_pairs[4]
    assert amount_pair_seq5 == 35
    assert '(()((.))(()(((.)((((.(().).)()))())((()(().).)(((..))(()))))((((((())))))())))).' in structures_seq5

    structures_seq6 = structures[5]
    amount_pair_seq6 = amount_pairs[5]
    assert amount_pair_seq6 == 63
    assert '(((()((((((.((()((((((()(((((())(().((())..))).)()((((...).))(()())(())())((.)))()()())))))(.)))))))(()((()(((((.)))))).))(((..))())).)))()..)..))' in structures_seq6



    