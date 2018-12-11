from prakt.scoring_func_parser.scoring_func_parser import ScoringMatrix


def test_scoring_matrix():
    scoring_matrix = ScoringMatrix("data/blosum62.txt", False, 2)
    # scoring_matrix.print_scoring_matrix()

    # test anything
    assert scoring_matrix.score("A", "A") == 4

    # test symmetric
    assert scoring_matrix.score("V", "N") == -3
    assert scoring_matrix.score("N", "V") == -3

    # gap score
    assert scoring_matrix.score("*", "A") == -4
    