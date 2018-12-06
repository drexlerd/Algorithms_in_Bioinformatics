from xpgma import XPGMA, Node

def test_example():
    m = [[0, 3, 12, 12, 9, 0, 0, 0, 0],
         [0, 0, 13, 13, 10, 0, 0, 0, 0],
         [0, 0, 0, 6, 7, 0, 0, 0, 0],
         [0, 0, 0, 0, 7, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0]
        ]
    l = [0, 1, 2, 3, 4]
    n = dict([(0, Node()), (1, Node()), (2, Node()), (3, Node()), (4, Node())])

    xpgma = XPGMA()

    root_node, node_dict = xpgma.generate_upgma(m, l, n)
    print(node_dict)