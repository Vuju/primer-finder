import numpy as np


def smith_waterman(primer, read, gap=-2, substitution_function=(lambda p, r: 2 if (p == r) else -1)):
    """
    An implementation of the Smith-Waterman algorithm for local sequence alignment.
    Extended with the functionality to jump triplets.
    :param primer: The first sequence. In my case a primer, hence the name.
    :param read: The second sequence. In this case the read in which I want to find my primer sequence.
    :param gap: The penalty value for skipping a triplet. Default -2.
    :param substitution_function: The reward and penalty values from comparing two characters.
    Has to be a function that accepts two characters to compare and returns a numerical value.
    By default, it returns +2 for a match, -1 for a mismatch.
    :return max_score: the highest score in the final score-matrix
    :return aligned_primer:the primer sequence with potential skips
    :return aligned_read: the read sequence with potential skips
    :return alignment_start_index_in_read: the starting index of the alignment within the read

    """
    # construction
    rows = len(primer) + 3
    cols = len(read) + 3

    score_matrix = np.zeros((rows, cols))
    traceback_matrix = np.zeros((rows, cols))

    max_score = 0
    max_pos = (0, 0)

    for i in range(3, rows):
        current_row = score_matrix[i, :]
        prev_row = score_matrix[i - 1, :]
        prev_prev_row = score_matrix[i - 3, :]

        for j in range(3, cols):
            match_score = prev_row[j - 1] + substitution_function(primer[i - 3], read[j - 3])
            delete = prev_row[j] + (2 * gap)
            insert = current_row[j - 1] + (2 * gap)

            # calculate triplet jump value
            del3 = prev_prev_row[j] + gap
            ins3 = current_row[j - 3] + gap

            scores = [0, match_score, delete, insert, del3, ins3]
            directions = np.argmax(scores)
            score_matrix[i][j] = scores[directions]
            traceback_matrix[i][j] = directions

            if scores[directions] > max_score:
                max_score = scores[directions]
                max_pos = (i, j)
        score_matrix[i, :] = current_row

    # traceback
    aligned_primer = []
    aligned_read = []
    i, j = max_pos

    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        match traceback_matrix[i, j]:
            case 1:  # Match/Mismatch
                aligned_primer.append(primer[i - 3])
                aligned_read.append(read[j - 3])
                i -= 1
                j -= 1
            case 2:  # Deletion
                aligned_primer.append(primer[i - 3])
                aligned_read.append('-')
                i -= 1
            case 3:  # Insertion
                aligned_primer.append('-')
                aligned_read.append(read[j - 3])
                j -= 1
            case 4:  # Triplet Deletion
                aligned_primer.append(primer[i - 5])
                aligned_read.append('---')
                i -= 3
            case 5:  # Triplet Insertion
                aligned_primer.append('---')
                aligned_read.append(read[j - 5])
                j -= 3
            case default:
                raise Exception(f"Smith-Waterman traceback failed! ({default})")

    aligned_primer = ''.join(reversed(aligned_primer))
    aligned_read = ''.join(reversed(aligned_read))
    alignment_start_index_in_read = j - 3

    return max_score, aligned_primer, aligned_read, alignment_start_index_in_read
