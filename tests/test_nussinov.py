from nussinov import Nussinov
from prakt.basis_classes.cell import Cell

def test_fill_matrix():
    nussinov = Nussinov()

    sequence = "GGUCCAC"

    d = [[Cell(i, j) for j in range(len(sequence) + 1)] for i in range(len(sequence))]

    nussinov.fill_matrix(d, sequence, loop_length=1)


def test_compute_optimal_abstract_structure():
    nussinov = Nussinov()

    sequence = "GGUCCAC"

    abstract_structures, amount = nussinov.compute_optimal_abstract_structure(sequence, complete_traceback=True)

    assert amount == 2

    assert abstract_structures == [(0, None), (0, (5, 1)), (1, (4, 2)), (0, (7, 1)), (3, (4, 2)), (3, (5, 2)), (3, (6, 3)), (0, (7, 2)), (7, (6, 3))]
    # e.g. (0, None), (0, (5, 1)), (1, (4, 2)) => ((.))..


def test_convert_abstract_structure_to_structure():
    sequence = "GGUCCAC"

    abstract_structures = [(0, None), (0, (5, 1)), (1, (4, 2)), (0, (7, 1)), (3, (4, 2)), (3, (5, 2)), (3, (6, 3)), (0, (7, 2)), (7, (6, 3))]

    nussinov = Nussinov()

    structures = nussinov.convert_abstract_structure_to_structure(sequence, 2, abstract_structures)

    assert structures == ['.((..))', '(.(..))', '((..).)', '((.)..)', '((.))..']


def test_nussinov():
    nussinov = Nussinov()

    structures, amount_pairs = nussinov.run("data/nussinov_guideline.fasta", True)

    print(structures)

    # given: (((((()())(.((()..))..)...))))
    assert structures ==  [['(.((((((.))(((((.)))))))(.))))', '(((.)(((.))(((((.))))))((.))))', '(((.)(((.)((((((.)))))))(.))))', '(((.)(((.))(((((.)))((.)))))))', '((.(((((.))(((((.)))))))(.))))', '(((.((((.))(((((.)))))))(.))))', '((((((.).)((((((.))))))(.)))))', '((((.(((.))(((((.))))))(.)))))', '(((((.((.)))((((.)))((.)))))))', '((((((.).)((((((.)))))(.))))))', '((((((.))(.(((((.)))))(.))))))', '((((((.))((.((((.)))))(.))))))', '((((((.))(((.(((.)))))(.))))))', '((((((.))(((((.)(.))))(.))))))', '((((((.))((((.((.)))))(.))))))', '((((((.))(((((.(.)))))(.))))))', '((((((.))(((((((.))))))..)))))', '((((.(((.))(((((.)))))(.))))))', '((((.(((..((((((.)))))))))))))', '((((.(((.)((((((.)))))).))))))', '((((.(((.)((((((.))))))).)))))', '(((((.((.)((((((.)))))))).))))', '((((((((.))(((((.)))))).).))))', '((((((((.))(((((.)))))))..))))', '(((((.((.)((((((.))))))))).)))', '((((((((.))(((((.)))))).)).)))', '((((((((.))(((((.))))))).).)))', '(((((.((.)((((((.)))))))))).))', '((((((((.))(((((.)))))).))).))', '((((((((.))(((((.))))))).)).))', '(((((.((.)((((((.))))))))))).)', '((((((((.))(((((.)))))).)))).)', '((((((((.))(((((.))))))).))).)']]

    assert '(((((()())(.((()..))..)...))))' in structures is True
