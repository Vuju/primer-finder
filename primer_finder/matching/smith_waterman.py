from primer_finder.config.config_loader import get_config_loader
from primer_finder.matching.dtos.match_result_dto import MatchResultDTO


def default_substitution_function(letter_in_primer, letter_in_read) -> float:
    """
    The default substitution function will return 2 if the primer character could be substituted by the sequence character match, -1 otherwise.

    :param letter_in_primer: A character from a primer sequence.
    :param letter_in_read: A character from a DNA sequence.

    :return: 2 for a match, -1 otherwise.
    """
    match letter_in_read:
        case 'A':
            return 2 if letter_in_primer in "AWMRDHVN" else -1
        case 'C':
            return 2 if letter_in_primer in "CSMYBHVN" else -1
        case 'G':
            return 2 if letter_in_primer in "GSKRBDVN" else -1
        case 'T':
            return 2 if letter_in_primer in "TWKYBDHN" else -1
        case _:
            raise Exception(f"unknown literal in read sequence: '{letter_in_read}'")


class SmithWaterman:
    """
    This class instantiates an implementation of the Smith-Waterman algorithm.
    Major changes are:
     - There is a separate penalty for triplet gaps, reasoning being that a triplet gap doesn't shift the reading frame.
     - There is an optional score bonus at the end of the dna sequence to accommodate for partial matches.
    """

    def __init__(self,
        gap_penalty = None,
        triplet_gap_penalty = None,
        end_of_read_bonus_value = None,
        sub_function = None
    ):
        """
        A configured SmithWaterman instance. As many parameters will be kept
        the same over many calls, they are to be set at construction time.

        :param gap_penalty: The penalty value for skipping a single nucleotide.
        :param triplet_gap_penalty: The penalty value for skipping a triplet.
        :param end_of_read_bonus_value: Value of bonus per missing base pair at the end of a read.
        :param sub_function: The reward and penalty values from comparing two characters.
        Has to be a function that accepts two characters to compare and returns a numerical value.
        Set to None for a default function.
        """
        # todo explore parameter space
        config = get_config_loader().get_config()

        self.gap_penalty = gap_penalty if gap_penalty is not None else config["algorithm"]["gap_penalty"]
        self.triplet_gap_penalty = triplet_gap_penalty if triplet_gap_penalty is not None else config["algorithm"]["triplet_gap_penalty"]
        self.end_of_read_bonus_value = end_of_read_bonus_value if end_of_read_bonus_value is not None else config["algorithm"]["end_of_read_bonus"]

        self.substitution_function = sub_function or default_substitution_function
        self.match_value = self.substitution_function("A", "A")

    def align_partial(self, primer, super_sequence, search_interval):
        """
        Align the primer sequence to the super sequence. You can set a search interval within the super sequence.
        :param primer: A primer sequence.
        :param super_sequence: A (longer) DNA sequence.
        :param search_interval: The search interval within the super sequence.
        :return:
        """
        match = self.align(
            primer_sequence=primer,
            dna_sequence=super_sequence[search_interval[0]:search_interval[1]],
            are_ends_eligible_for_bonus=(search_interval[0] == 0, search_interval[1] == len(super_sequence))
        )
        match.start_index += search_interval[0]
        match.end_index += search_interval[0]
        return match


    def align(
            self,
            primer_sequence,
            dna_sequence,
            are_ends_eligible_for_bonus
    ) -> MatchResultDTO:
        """
        An implementation of the Smith-Waterman algorithm for local sequence alignment.
        Extended with the functionality to jump triplets.

        :param primer_sequence: The first sequence. In my case a primer, hence the name.
        :param dna_sequence: The second sequence. In this case the read in which I want to find my primer sequence.
        :param are_ends_eligible_for_bonus: Tuple whether to give a bonus for a partial match at the beginning and end of read.

        :returns: A MatchResult object containing the match score, traceback within the dna_sequence
        as well as indexes of the match start and end within the dna_sequence.
        """
        # construction
        rows = len(primer_sequence) + 3
        cols = len(dna_sequence) + 3

        score_matrix = [[0] * cols for _ in range(rows)]
        # For traceback: 0=none, 1=diagonal(match), 2=up(deletion), 3=left(insertion), 4=up3(triplet deletion), 5=left3(triplet insertion)
        traceback_matrix = [[0] * cols for _ in range(rows)]

        max_score = 0
        max_pos = (0, 0)

        # for beginning (and end later), define a bonus to encourage matching partial primers
        if are_ends_eligible_for_bonus[0]:
            for i in range(2, rows):
                for j in range(0, 3):
                    score_matrix[i][j] = self.end_of_read_bonus_value * (i - 2)

        for i in range(3, rows):
            for j in range(3, cols):
                # Calculate all possible scores
                match_score = score_matrix[i - 1][j - 1] + (self.substitution_function(primer_sequence[i - 3], dna_sequence[j - 3]))
                delete = score_matrix[i - 1][j] + self.gap_penalty
                insert = score_matrix[i][j - 1] + self.gap_penalty
                del3 = (score_matrix[i - 3][j] + self.triplet_gap_penalty)
                ins3 = (score_matrix[i][j - 3] + self.triplet_gap_penalty)

                # Find maximum score and corresponding direction
                scores = [0, match_score, delete, insert, del3, ins3]
                max_index = scores.index(max(scores))

                score_matrix[i][j] += scores[max_index]
                traceback_matrix[i][j] = max_index if scores[max_index] > 0 else 0

                if score_matrix[i][j] > max_score:
                    max_score = score_matrix[i][j]
                    max_pos = (i, j)

        # for (beginning earlier and) end, define a bonus to encourage matching partial primers
        if are_ends_eligible_for_bonus[1]:
            last_column = cols - 1
            for i in range(3, rows):
                score_matrix[i][last_column] += max(0, self.end_of_read_bonus_value * (rows - i - 1))
                if score_matrix[i][last_column] > max_score:
                    max_score = score_matrix[i][last_column]
                    max_pos = (i, last_column)

        # perform traceback
        aligned_read = []
        i, j = max_pos


        # Continue until we hit a cell with zero score or reach the matrix boundary
        while i >= 3 and j >= 3 and score_matrix[i][j] > 0:
            direction = traceback_matrix[i][j]

            if direction == 0:  # End of alignment
                break

            elif direction == 1:  # Diagonal (match/mismatch)
                aligned_read.append(dna_sequence[j - 3])
                i -= 1
                j -= 1

            elif direction == 2:  # Up (deletion in read)
                aligned_read.append('-')
                i -= 1

            elif direction == 3:  # Left (insertion in read)
                aligned_read.append(dna_sequence[j - 3])
                j -= 1

            elif direction == 4:  # Up 3 (triplet deletion)
                for k in range(3):
                    if i - k > 2:
                        aligned_read.append('-')
                i -= 3

            elif direction == 5:  # Left 3 (triplet insertion)
                for k in range(3):
                    if j - k > 2:
                        aligned_read.append(dna_sequence[j - 3 - k])
                j -= 3

        aligned_read = ''.join(reversed(aligned_read))

        # There were 5 in 2.3M sequences, which managed to set their index to -1 instead of 0.
        # Since their alignment worked otherwise, I set the index to min 0. Maybe their score is also off by 1 but who cares.
        alignment_start_index_in_read = max(0, j - 2)

        return MatchResultDTO(max_score, aligned_read, alignment_start_index_in_read, alignment_start_index_in_read + len(aligned_read))
