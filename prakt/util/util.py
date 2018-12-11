import math


def Min(sequence):
    """
    Args:
      sequence (list(tuples)): A list of tuples where the first item of a tuple is the value of the item
                               and second item is the case identifier
    
    Returns:
      list(tuple): A list of tuples with minimum value and the corresponding case identifier
    """
    min_value = min(sequence, key=lambda x: x[0])[0]
    return [t for t in sequence if t[0] == min_value]


def Max(sequence):
    """
    Args:
      sequence (list(tuples)): A list of tuples where the first item of a tuple is the value of the item
                               and second item is the case identifier
    
    Returns:
      tuple: A list of tuples with maximum value and the corresponding case identifier
    """
    max_value = max(sequence, key=lambda x: x[0])[0]
    return [t for t in sequence if t[0] == max_value]


def compute_traceback(start, all=False):
    """Compute the alignment(s) generated by tracing back in the dp matrix

    Args:
        start (cell): the cell from where to start dfs
        all (boolean): Indicates if all optimal alignments should get returned
                        or just one
    
    Return:
        list(list(Cell)) : A list containing all paths of integers which yield an optimal alignment
                           The first element is the bottom right case.
                           The last element is the 1 1 case.
    """
    tracebacks = compute_traceback_dfs(start)
    if all:
        return tracebacks
    return [tracebacks[0]]


def compute_traceback_dfs(current_cell, current_path=[]):
    """Recursively creates the list of all tracebacks.

    Initial call: Bottom right cell in the dynamic programming matrix

    Args:
        current_cell (Cell) : The cell where we are currently in the traceback
        current_path (list(int)) : List of cases that are on the traceback path
                                   first element is case in bottom right cell
                                   last element is the previous cases
    """
    # base case: leaf
    if not current_cell.pre:
        return [current_path]
    # recursive case: inner node
    tracebacks = []
    for pre in current_cell.pre:
        # current_path.append(pre)
        # Note: current_path + [pre] instead of append because he have to copy the path
        tracebacks.extend(compute_traceback_dfs(pre[0], current_path + [pre]))
    return tracebacks


def similarity_to_distance(nw, pairwise_alignment, scoring_matrix):
    """Converts similarity score from a pairwise alignment to a distance score
    using approximation algorithm
    
    D(a,b) = - log(S_{a,b}^{eff})

    S_{a,b}^{eff} = (S(a,b) - S_{rand}) / (S_{a,b}^{max} - S_{rand})

    S_{rand} = (1/|A|) * (sum_{x,y in Sigma x Sigma} S(x,y) * N_a(x) * N_b(y)) + gaps(A) * S(-,*)

    S_{a,b}^{max} = (S(a,a) + S(b,b)) / 2
    """
    seq1 = pairwise_alignment[0].replace("_", "")
    seq2 = pairwise_alignment[1].replace("_", "")

    S_ab, _ = nw.compute_optimal_alignments(seq1, seq2, scoring_matrix, complete_traceback=False)

    S_aa, _ = nw.compute_optimal_alignments(seq1, seq1, scoring_matrix, complete_traceback=False)

    S_bb, _ = nw.compute_optimal_alignments(seq2, seq2, scoring_matrix, complete_traceback=False)

    S_ab_max = (S_aa + S_bb) / 2

    S_rand = (1 / len(pairwise_alignment[0])) * \
        sum([scoring_matrix.score(scoring_matrix.alphabet[i], scoring_matrix.alphabet[j]) * count_occurences_symbol_in_seq(seq1, scoring_matrix.alphabet[i]) * count_occurences_symbol_in_seq(seq2, scoring_matrix.alphabet[j]) for i in range(len(scoring_matrix.alphabet)) for j in range(len(scoring_matrix.alphabet))]) + count_gaps_in_pairwise_alignment(pairwise_alignment) * scoring_matrix.cost_gap_open

    S_eff = (S_ab - S_rand) / (S_ab_max - S_rand)

    #print(pairwise_alignment)
    #print("S_ab %5.2f" % S_ab)
    #print("S_ab_max %5.2f" % S_ab_max)
    #print("S_rand %5.2f" % S_rand)
    #print("S_eff %5.2f" %S_eff)

    return - math.log(S_eff)


def count_occurences_symbol_in_seq(seq, symbol):
    """Counts the number of occurences of symbol in seq

    Args:
        seq (str): A sequence
        symbol (char): A character
    """
    count = 0
    for x in seq:
        if x == symbol:
            count += 1
    return count


def count_gaps_in_pairwise_alignment(pairwise_alignment):
    """Counts the number of gaps in the given pairwise alignment. Since Feng-Dootlittle
    exchanges gaps with special character "X" we also need to take care about this special
    character

    Args:
        pairwise_alignment (list(str)): A 2 element list containing a gapped sequences.
                                        Therefore the lengths are equal
                                        
    """
    count = 0
    for i in range(len(pairwise_alignment[0])):
        if pairwise_alignment[0][i] in ["_", "X"] or pairwise_alignment[1][i] in ["_", "X"]:
            count += 1
    return count
