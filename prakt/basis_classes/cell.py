from enum import Enum


class Cell(object):
    """This class represents a cell in the dynamic programming matrix
    """
    def __init__(self, i, j):
        # position
        self.i = i
        self.j = j
        # list of tuples containing predecessor with the respective case how it was reached
        # used in traceback for 
        self.pre = []
        # the value
        self.value = 0

    def AddPredecessor(self, pre_cell, case):
        self.pre.append((pre_cell, case))

    def SetValue(self, value):
        self.value = value

    def GetValue(self):
        return self.value

    def get_position(self):
        return self.i, self.j


class MatrixType(Enum):
    D = 0
    P = 1
    Q = 2