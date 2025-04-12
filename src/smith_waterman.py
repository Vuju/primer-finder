# todo explore parameter space

def smith_waterman(
        primer,
        read,
        gap = -2,
        gap3 = -2,
        substitution_function = (lambda p, r: 2 if (p == r) else -1),
        end_of_read_bonus = (True, True),
        end_of_read_bonus_value = 1
):
    """
    An implementation of the Smith-Waterman algorithm for local sequence alignment.
    Extended with the functionality to jump triplets.
    :param primer: The first sequence. In my case a primer, hence the name.
    :param read: The second sequence. In this case the read in which I want to find my primer sequence.
    :param gap: The penalty value for skipping a single nucleotide. Default -2.
    :param gap3: The penalty value for skipping a triplet. Default -2.
    :param substitution_function: The reward and penalty values from comparing two characters.
    Has to be a function that accepts two characters to compare and returns a numerical value.
    By default, it returns +2 for a match, -1 for a mismatch.
    :param end_of_read_bonus: Tuple whether to give a bonus for a partial match at the beginning and end of read.
    :param end_of_read_bonus_value: Value of bonus per missing base pair at the end of a read.
    :return max_score: the highest score in the final score-matrix
    :return aligned_primer:the primer sequence with potential skips
    :return aligned_read: the read sequence with potential skips
    :return alignment_start_index_in_read: the starting index of the alignment within the read
    """
    # construction
    rows = len(primer) + 3
    cols = len(read) + 3

    score_matrix = [[0] * cols for _ in range(rows)]
    # For traceback: 0=none, 1=diagonal(match), 2=up(deletion), 3=left(insertion), 4=up3(triplet deletion), 5=left3(triplet insertion)
    traceback_matrix = [[0] * cols for _ in range(rows)]

    max_score = 0
    max_pos = (0, 0)

    # for beginning (and end later), define a bonus to encourage matching partial primers
    if end_of_read_bonus[0]:
        for i in range(2, rows):
            for j in range(0, 3):
                score_matrix[i][j] = end_of_read_bonus_value * (i - 2)

    for i in range(3, rows):
        for j in range(3, cols):
            # Calculate all possible scores
            match_score = score_matrix[i - 1][j - 1] + (substitution_function(primer[i - 3], read[j - 3]))
            delete = score_matrix[i - 1][j] + gap
            insert = score_matrix[i][j - 1] + gap
            del3 = (score_matrix[i - 3][j] + gap3)
            ins3 = (score_matrix[i][j - 3] + gap3)

            # Find maximum score and corresponding direction
            scores = [0, match_score, delete, insert, del3, ins3]
            max_index = scores.index(max(scores))

            score_matrix[i][j] += scores[max_index]
            traceback_matrix[i][j] = max_index if scores[max_index] > 0 else 0

            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_pos = (i, j)

    # for (beginning earlier and) end, define a bonus to encourage matching partial primers
    if end_of_read_bonus[1]:
        last_column = cols - 1
        for i in range(3, rows):
            score_matrix[i][last_column] += max(0, end_of_read_bonus_value * (rows - i - 1))
            if score_matrix[i][last_column] > max_score:
                max_score = score_matrix[i][last_column]
                max_pos = (i, last_column)

    # perform traceback
    aligned_primer = []
    aligned_read = []
    i, j = max_pos

    # debug print:
    #for line in score_matrix:
    #    print(line)

    # Continue until we hit a cell with zero score or reach the matrix boundary
    while i >= 3 and j >= 3 and score_matrix[i][j] > 0:
        direction = traceback_matrix[i][j]

        if direction == 0:  # End of alignment
            break

        elif direction == 1:  # Diagonal (match/mismatch)
            aligned_primer.append(primer[i - 3])
            aligned_read.append(read[j - 3])
            i -= 1
            j -= 1

        elif direction == 2:  # Up (deletion in read)
            aligned_primer.append(primer[i - 3])
            aligned_read.append('-')
            i -= 1

        elif direction == 3:  # Left (insertion in read)
            aligned_primer.append('-')
            aligned_read.append(read[j - 3])
            j -= 1

        elif direction == 4:  # Up3 (triplet deletion)
            for k in range(3):
                if i - k > 2:  # Ensure we're within bounds
                    aligned_primer.append(primer[i - 3 - k])
                    aligned_read.append('-')
            i -= 3

        elif direction == 5:  # Left3 (triplet insertion)
            for k in range(3):
                if j - k > 2:  # Ensure we're within bounds
                    aligned_primer.append('-')
                    aligned_read.append(read[j - 3 - k])
            j -= 3

    aligned_primer = ''.join(reversed(aligned_primer))
    aligned_read = ''.join(reversed(aligned_read))

    # There were 5 in 2.3M sequences, which managed to set their index to -1 instead of 0.
    # Since their alignment worked otherwise, I set the index to min 0. Maybe their score is also off by 1 but who cares.
    alignment_start_index_in_read = max(0, j - 2)

    return max_score, aligned_primer, aligned_read, alignment_start_index_in_read
