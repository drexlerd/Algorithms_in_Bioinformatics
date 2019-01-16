"""
Module name: Nussinov
Module author: dominik drexler <drexlerd@informatik.uni-freiburg.de>
"""

from prakt.nussinov import NussinovBase
from prakt.fasta_parser.fasta_parser import parse_fasta, check_sequences_alphabet, SequenceType
from prakt.basis_classes.cell import Cell
from enum import Enum
from prakt.util.util import Max, compute_traceback
import argparse

class Case(Enum):
    UNPAIRED = 0,
    PAIRED_SUCCESS = 1,
    PAIRED_FAIL = 2


@NussinovBase.register
class Nussinov(NussinovBase):
    def _is_base_pair(self, c1, c2):
        """Returns True if c1 and c2 build a base pair and False otheriwse

        Args:
          c1 (char): A character
          c2 (char): A character
        """
        if c1 == "A" and c2 == "U" or c1 == "U" and c2 == "A" \
          or c1 == "G" and c2 == "C" or c1 == "C" and c2 == "G" \
          or c1 == "G" and c2 == "U" or c1 == "U" and c2 == "G":
            return True
        return False

    
    def fill_matrix(self, d, sequence, loop_length=1):
        """Inductive case of DP algorithm.

        Args:
          d (list(list(Cell))): The DP matrixs initialized to zero values
          sequence (str): The given string
        """
        min_loop_length = 2
        for l in range(min_loop_length, len(sequence) + 1):
            for i in range(1, len(sequence) - l + 1):
                j = i + l
                # print("%d, %d" % (i,j))
                # i and j are the indices iterated in diagonal
                case_dist = [(d[i][j-1].value, (i, j, None), Case.UNPAIRED)]
                
                for k in range(i, j - loop_length):
                    if self._is_base_pair(sequence[k-1],sequence[j-1]):
                        # print("%d, %d, %s, %s" % (i, j, sequence[k-1], sequence[j-1]))
                        case_dist.append((d[i][k-1].value + d[k+1][j-1].value + 1, (i, j, k), Case.PAIRED_SUCCESS))
                    else:
                        case_dist.append((0, (i, j, k), Case.PAIRED_FAIL))

                result_sequence = Max(case_dist)

                # store result of the recursion
                d[i][j].SetValue(result_sequence[0][0])  # set value
                # set predecessors
                for result_case in result_sequence:
                    d[i][j].AddPredecessor(result_case[1], result_case[2])
        
        #print()
        #for i in range(len(sequence)):
        #    for j in range(len(sequence) + 1):
        #        print("%6.2f" % (d[i][j].value), end='')
        #    print()


    def traceback_rec(self, d, i, j, structure_index=0, abstract_structures=[(0, None)]):
        """Computes abstract treelike structure representation of all optimal pairing.
                
        Initial call: traceback_rec(d, 0, len(sequence), 0, [(0, None)])

        Args:
          d (list(list(Cell))): DP matrix (with traceback pointers)
          i (int): smallest index in sequence (included)
          j (int): biggest index in sequence (included)
          structure_index (int): the index in abstract_structures we are currently at in our recursion
          abstract_structures (list(tuple(int,tuple(int,int)))): pointer to head recursion + a tuples of base pairings

        Returns:
          abstract_structures (list(tuple(int,tuple(int,int)))): pointer to head recursion + a tuples of base pairings
        """
        current_cell = d[i][j]
        if not self._is_traceback_base_case(current_cell):
            for (i, j, k), case in current_cell.pre:
                # copy the structure because of multiple predecessor
                if case == Case.UNPAIRED:
                    self.traceback_rec(d, i, j-1, structure_index, abstract_structures)
                elif case == Case.PAIRED_SUCCESS:
                    abstract_structures.append((structure_index, (j, k)))
                    self.traceback_rec(d, i, k-1, len(abstract_structures) - 1, abstract_structures)
                    self.traceback_rec(d, k + 1, j - 1, len(abstract_structures) - 1, abstract_structures)
                elif case == Case.PAIRED_FAIL:
                    self.traceback_rec(d, i, k-1, len(abstract_structures) - 1, abstract_structures)
                    self.traceback_rec(d, k + 1, j - 1, len(abstract_structures) - 1, abstract_structures)
            return abstract_structures


    def _is_traceback_base_case(self, cell):
        """Helper function of traceback_rec.
        Returns True if the cell of the initialization was reached

        Args:
          cell (Cell): A cell for which to test if its a base case cell
        """
        if cell.i + 2 > cell.j:
            return True
        return False


    def compute_optimal_abstract_structure(self, sequence, complete_traceback):
        """For a given sequence compute the abstract structure

        An abstract structure is equivalent to the resulting structure,
        but the representation differs because of the way its computed

        Args:
          sequence (str): A string
          complete_traceback (bool): True, full traceback

        Returns:
          abstract_structures (list(tuple(int,tuple(int,int)))): pointer to head recursion + a tuples of base pairings
          amount_pairs (int): The amount of base pair in optimal abstract structure
        """
        d = [[Cell(i, j) for j in range(len(sequence) + 1)] for i in range(len(sequence) + 1)]

        self.fill_matrix(d, sequence, loop_length=1)

        abstract_structures = self.traceback_rec(d, 1, len(sequence), 0, [(0, None)]) 

        return abstract_structures, d[1][len(sequence)].value
        

    def convert_abstract_structure_to_structure(self, sequence, amount_pairs, abstract_structures):
        """Converts an abstract treelike structure into a human readable form.

        Args:
          sequence (str): A string
          amount_pairs (int): The amount of base pair in optimal abstract structure
          abstract_structures (list(tuple(int,tuple(int,int)))): pointer to head recursion + a tuples of base pairings
        """
        structures = []
        # track which abstract structures are already consumed
        for index, (j, k) in reversed(abstract_structures[1:]):
            count = 0
            r = ["."] * len(sequence)
            while True:
                r[j-1] = ")"
                r[k-1] = "("
                count += 1
                if count == amount_pairs:
                    structures.append("".join(r))
                if abstract_structures[index][1] == None:
                    break
                index, (j, k) = abstract_structures[index]
        return structures


    def run(self,
        seq_fasta_fn,
        complete_traceback):
        """Given a fasta file computes optimal structures using the nussinov algorithm
        If complete traceback then all optimal structures are returned, else 1

        Args:
          seq_fast_fn (str): A fasta file
          complete_traceback (bool): All optimal structures, if True, else 1
        """
        # sequences with their ids
        records = parse_fasta(seq_fasta_fn)

        results = []
  
        for r in records:
            sequence = str(r.seq)

            abstract_structures, amount_pairs = self.compute_optimal_abstract_structure(sequence, complete_traceback)

            structures = self.convert_abstract_structure_to_structure(sequence, amount_pairs, abstract_structures)

            results.append(structures)

        return results, amount_pairs



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Nussinov command line tool")
    parser.add_argument("seq_fasta_fn", type=str)
    parser.add_argument("--c", "--complete_traceback", action='store_true')
    args = parser.parse_args()

    nussinov = Nussinov()

    # sequences with their ids
    records = parse_fasta(args.seq_fasta_fn)
    complete_traceback = args.c

    results, amount_pairs = nussinov.run(args.seq_fasta_fn, complete_traceback)


    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print("Nussinov Results")
    print("Maximal number of base pairs: %d" % amount_pairs)
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    for i in range(len(results)):
        sequence = str(records[i].seq)
        seqID = records[i].seq

        structures = results[i]

        print("Optimal structures for sequence %s" % (seqID))
        print("Total optimal structures: %d" % (len(structures)))
        print(sequence)
        for structure in structures:
            print(structure)
        print()