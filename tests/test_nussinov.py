from nussinov import Nussinov
from prakt.basis_classes.cell import Cell

def test_fill_matrix():
    nussinov = Nussinov()

    sequence = "GGUCCAC"

    d = [[Cell(i, j) for j in range(len(sequence) + 1)] for i in range(len(sequence))]

    nussinov.fill_matrix(d, sequence, loop_length=1)


def test_compute_optimal_structure():
    nussinov = Nussinov()

    sequence = "GGUCCAC"

    structures = nussinov.compute_optimal_structure(sequence, complete_traceback=True)

    print(structures)