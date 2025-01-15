import numpy as np


min_loop_length = 4


def complementary(nucleotide1, nucleotide2):
    """ Return True if RNA nucleotides are Watson-Crick base pairs. """
    complement = {"A": "U", "U": "A", "G": "C", "C": "G"}
    if complement[nucleotide1] == nucleotide2:
        return True
    return False


def nussinov(rna):
    """ Calculate structure matrix as per the Nussinov algorithm. """
    n = len(rna)

    # dp[i][j] stores the score of the optimal pairing between indices i and j
    dp = np.zeros((n, n))

    # fill the dp matrix diagonally
    for k in range(n):
        for i in range(n - k):
            j = i + k

            # case 1: check if j not paired
            unpaired = dp[i][j-1]

            # case 2: check if j can be involved in a pairing with a position t
            pairing = [0]
            for t in range(i, j - min_loop_length):
                if complementary(rna[t], rna[j]):
                    pairing += [dp[i][t-1] + dp[t+1][j-1] + 1]
            paired = max(pairing)

            dp[i][j] = max(unpaired, paired)

    return dp


def traceback(i, j, dp, rna, structure):
    """ Find the optimal substructure between indices [i...j]. """
    if i >= j:
        return

    # case 1: j is unpaired
    if dp[i][j] == dp[i][j-1]:
        traceback(i, j-1, dp, rna, structure)
        return

    # case 2: pairing j with a matching index k to its left
    for k in range(i, j - min_loop_length):
        if complementary(rna[k], rna[j]):

            # handle boundary value
            if k == 0:
                if dp[i][j] == dp[k+1][j-1] + 1:
                    structure.append((k, j))
                    traceback(k + 1, j - 1, dp, rna, structure)
                    return

            if dp[i][j] == dp[i][k-1] + dp[k+1][j-1] + 1:
                # add the pair (k, j) to list of pairs
                structure.append((k, j))

                # move the recursion to the two substructures formed by this pairing
                traceback(i, k-1, dp, rna, structure)
                traceback(k+1, j-1, dp, rna, structure)
                return


def write_structure(rna, structure):
    """ Represent the structure in dot-bracket notation. """
    dot_bracket = ["." for _ in range(len(rna))]
    for s in structure:
        dot_bracket[min(s)] = "("
        dot_bracket[max(s)] = ")"
    return "".join(dot_bracket)


if __name__ == "__main__":
    S = "UAAGCGCAUAGACGAUGCAACCGCAUCGCCUAGCGCUUAGAUCGGCCAACAGCAUCGGUUGGCCAAUAGGCUGCACGCCAUCAUUUGCAUGGUGUGUAAAGCUU"
    N = nussinov(S)
    structure = []
    traceback(0, len(S)-1, N, S, structure)
    print(write_structure(S, structure))
