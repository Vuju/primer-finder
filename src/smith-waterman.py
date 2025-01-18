
# Plan 2: reformat input parameters
# Plan 3: rework the traceback, especially the output


def smith_waterman(primer, read, gap=-2, substitution_matrix=(lambda x, y: 2 if (x == y) else -1)):

    # construction
    rows = len(primer) + 1
    cols = len(read) + 1
    score_matrix = [[0] * cols for _ in range(rows)]

    max_score = 0
    max_pos = (0, 0)

    for i in range(1, rows):
        for j in range(1, cols):
            match_score = score_matrix[i - 1][j - 1] + (substitution_matrix(primer[i - 1], read[j - 1]))
            delete = score_matrix[i - 1][j] + (2 * gap)
            insert = score_matrix[i][j - 1] + (2 * gap)

            # calculate triplet jump value
            del3 = (score_matrix[i - 3][j] + gap) if (i - 3 >= 0) else 0
            ins3 = (score_matrix[i][j - 3] + gap) if (j - 3 >= 0) else 0

            score_matrix[i][j] = max(0, match_score, delete, insert, del3, ins3)

            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_pos = (i, j)

    # traceback from ChatGPT, Todo: rework
    aligned_primer = []
    aligned_read = []
    i, j = max_pos

    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        if score_matrix[i][j] == score_matrix[i - 1][j - 1] + (substitution_matrix(primer[i - 1], read[j - 1])):
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
        elif (i - 3 >= 0) & score_matrix[i][j] == score_matrix[i - 3][j] + gap:
            aligned_primer.append(primer[i - 3])
            aligned_read.append('---')
            i -= 3
        elif (j - 3 >= 0) & score_matrix[i][j] == score_matrix[i][j - 3] + gap:
            aligned_primer.append('---')
            aligned_read.append(read[j - 3])
            j -= 3
        else:
            raise Exception("Smith-Waterman traceback failed!")

    aligned_primer = ''.join(reversed(aligned_primer))
    aligned_read = ''.join(reversed(aligned_read))
    alignment_start_index_in_read = j

    return max_score, aligned_primer, aligned_read, alignment_start_index_in_read
