"""
author: dominik drexler <drexlerd@informatik.uni-freiburg.de>
"""
from enum import Enum


class MetricType(Enum):
    SIMILARITY = 0
    DISTANCE = 1


class ScoringMatrix(object):
    """This class represents a scoring function.

    It has a function score to compute the score for mutating from one symbol
    to another symbol.
    """
    def __init__(self, filename, is_distance_fn):
        # the quadratic scoring matrix to retrieve scores from.
        self.scoring_matrix = []
        # mapping of symbol from the alphabet to an index, row or column in the matrix
        # Note: dont assume same ordering of symbol in rows / cols
        self.symbol_to_position_col = dict()
        self.symbol_to_position_row = dict()
        self.position_to_symbol_col = dict()
        self.position_to_symbol_row = dict()

        # alphabet
        self.alphabet = []

        # the type of the scoring metric (distance or similarity)
        self.metric_type = MetricType.SIMILARITY

        self._parse_scoring_func(filename)

        self._set_metric_type(is_distance_fn)

        # self._add_neutral_symbol()  # added special case in nw for that

    
    def _add_neutral_symbol(self):
        self.symbol_to_position_col["X"] = len(self.symbol_to_position_col)
        self.symbol_to_position_row["X"] = len(self.symbol_to_position_row)
        for i in range(len(self.scoring_matrix[0])):
            self.scoring_matrix[i].append(0)
        self.scoring_matrix.append([0] * len(self.symbol_to_position_col))


    def _parse_scoring_func(self, filename):
        """Fill the symbol_to_position_col, symbol_to_position_row dictionaries and the scoring_matrix
        """
        path = filename
        with open(path, 'r') as file:
            i = 0  # count the number of lines of the matrix
            symbol_row_list = []
            symbol_col_list = []
            for line in file:
                if line.startswith("#"):  # this indicates a comment line
                    continue
                line_list = line.split()
                if i == 0:  # this is the first line containing only symbols
                    symbol_col_list = line_list
                    self.symbol_to_position_col = dict(zip(symbol_col_list, range(len(symbol_col_list))))
                    self.position_to_symbol_col = dict(zip(range(len(symbol_col_list)), symbol_col_list))
                    self.alphabet = symbol_col_list
                else:
                    line_list_scores = [int(i) for i in line_list[1:]]
                    line_symbol = line_list[0]
                    self.scoring_matrix.append(line_list_scores)  # first element is the symbol
                    symbol_row_list.append(line_symbol)
                i += 1
        self.symbol_to_position_row = dict(zip(symbol_row_list, range(len(symbol_row_list))))
        self.position_to_symbol_row = dict(zip(range(len(symbol_row_list)), symbol_row_list))
        assert symbol_col_list == symbol_row_list  # for now check for being quadratic

    
    def _set_metric_type(self, is_distance_fn):
        if is_distance_fn:
            self.metric_type = MetricType.DISTANCE
        else:
            self.metric_type = MetricType.SIMILARITY


    def score(self, symbol1, symbol2):
        """Returns the score of mutating from 'symbol1' to 'symbol2'
        """
        try:
            pos1 = self.symbol_to_position_row[symbol1]
            pos2 = self.symbol_to_position_col[symbol2]
        except KeyError:
            return None
        return self.scoring_matrix[pos1][pos2]


    def print_scoring_matrix(self):
        """Prints a the scoring matrix nicely
        """
        print(3 * " ", end='')
        for j in range(len(self.position_to_symbol_col)):
            print("%3s" % (self.position_to_symbol_col[j]), end='')
        print()
        for i in range(len(self.position_to_symbol_row)):
            print("%3s" % (self.position_to_symbol_row[i]), end='')
            for j in range(len(self.position_to_symbol_col)):
                print("%3s" % (self.scoring_matrix[i][j]), end='')
            print()
        print("end")

