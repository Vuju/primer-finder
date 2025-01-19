

def smith_waterman(primer, read, gap=-2, substitution_function=(lambda p, r: 2 if (p == r) else -1)):
    """
    An implementation of the Smith-Waterman algorithm for local sequence alignment.
    Extended with the functionality to jump triplets.
    :param primer: The first sequence. In my case a primer, hence the name.
    :param read: The second sequence. In this case the read in which I want to find my primer sequence.
    :param gap: The penalty value for skipping a triplet. Default -2.
    :param substitution_function: The reward and penalty values from comparing two characters.
    Has to be a function that accepts two characters to compare and returns a numerical value.
    By default it returns +2 for a match, -1 for a mismatch.
    :return max_score: the highest score in the final score-matrix
    :return aligned_primer:the primer sequence with potential skips
    :return aligned_read: the read sequence with potential skips
    :return alignment_start_index_in_read: the starting index of the alignment within the read

    """
    # construction
    rows = len(primer) + 1
    cols = len(read) + 1
    score_matrix = [[0] * cols for _ in range(rows)]

    max_score = 0
    max_pos = (0, 0)

    for i in range(1, rows):
        for j in range(1, cols):
            match_score = score_matrix[i - 1][j - 1] + (substitution_function(primer[i - 1], read[j - 1]))
            delete = score_matrix[i - 1][j] + (2 * gap)
            insert = score_matrix[i][j - 1] + (2 * gap)

            # calculate triplet jump value
            del3 = (score_matrix[i - 3][j] + gap) if (i - 3 >= 0) else 0
            ins3 = (score_matrix[i][j - 3] + gap) if (j - 3 >= 0) else 0

            score_matrix[i][j] = max(0, match_score, delete, insert, del3, ins3)

            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_pos = (i, j)

    # traceback
    aligned_primer = []
    aligned_read = []
    i, j = max_pos

    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        if score_matrix[i][j] == score_matrix[i - 1][j - 1] + (substitution_function(primer[i - 1], read[j - 1])):
            aligned_primer.append(primer[i - 1])
            aligned_read.append(read[j - 1])
            i -= 1
            j -= 1
        elif score_matrix[i][j] == score_matrix[i - 1][j] + (2 * gap):
            aligned_primer.append(primer[i - 1])
            aligned_read.append('-')
            i -= 1
        elif score_matrix[i][j] == score_matrix[i][j - 1] + (2 * gap):
            aligned_primer.append('-')
            aligned_read.append(read[j - 1])
            j -= 1
        elif (i - 3 >= 0) and (score_matrix[i][j] == score_matrix[i - 3][j] + gap):
            aligned_primer.append(primer[i - 3])
            aligned_read.append('---')
            i -= 3
        elif (j - 3 >= 0) and (score_matrix[i][j] == score_matrix[i][j - 3] + gap):
            aligned_primer.append('---')
            aligned_read.append(read[j - 3])
            j -= 3
        else:
            raise Exception("Smith-Waterman traceback failed!")

    aligned_primer = ''.join(reversed(aligned_primer))
    aligned_read = ''.join(reversed(aligned_read))
    alignment_start_index_in_read = j

    return max_score, aligned_primer, aligned_read, alignment_start_index_in_read
