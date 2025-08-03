import logging

from primer_finder.config.config_loader import get_config_loader
from primer_finder.matching.dtos.match_result_dto import MatchResultDTO

logger = logging.getLogger(__name__)

def default_substitution_function(letter_in_primer: str, letter_in_read: str) -> float:
    """
    The default substitution function will return 2 if the primer character could be substituted by the sequence character match, -1 otherwise.

    :param letter_in_primer: A character from a primer sequence.
    :param letter_in_read: A character from a DNA sequence.

    :return: 2 for a match, -1 otherwise.
    """
    if letter_in_primer is None or letter_in_read is None:
        logger.error(f"Invalid input to substitution function: primer={letter_in_primer}, read={letter_in_read}")
        return -1
    
    try:
        match letter_in_primer:
            case 'A':
                return 2 if letter_in_read in "A" else -1
            case 'C':
                return 2 if letter_in_read in "C" else -1
            case 'G':
                return 2 if letter_in_read in "G" else -1
            case 'T':
                return 2 if letter_in_read in "TU" else -1
            case 'U':
                return 2 if letter_in_read in "TU" else -1
            case 'W':
                return 2 if letter_in_read in "WATU" else -1
            case 'S':
                return 2 if letter_in_read in "SCG" else -1
            case 'M':
                return 2 if letter_in_read in "MAC" else -1
            case 'K':
                return 2 if letter_in_read in "KGTU" else -1
            case 'R':
                return 2 if letter_in_read in "RAG" else -1
            case 'Y':
                return 2 if letter_in_read in "YCTU" else -1
            case 'B':
                return 2 if letter_in_read in "BCGTSKYU" else -1
            case 'D':
                return 2 if letter_in_read in "DAGTWKRU" else -1
            case 'H':
                return 2 if letter_in_read in "HACTWMYU" else -1
            case 'V':
                return 2 if letter_in_read in "VACGSMR" else -1
            case 'N':
                return 2
            case '-':
                logger.warning("Gap character '-' found in sequence - this should never happen in eyeBOLD")
                return 0
            case _:
                logger.error(f"Unknown character in primer sequence: '{letter_in_primer}'")
                return -1
    except Exception as e:
        logger.error(f"Error in substitution function: {str(e)}")
        return -1


class SmithWaterman:
    """
    This class instantiates an implementation of the Smith-Waterman algorithm.
    Major changes are:
     - There is a separate penalty for triplet gaps, reasoning being that a triplet gap doesn't shift the reading frame.
     - There is an optional score bonus at the end of the dna sequence to accommodate for partial matches.
    """

    def __init__(self,
        gap_penalty: int = None,
        triplet_gap_penalty: int = None,
        end_of_read_bonus_value: int = None,
        sub_function: callable = None,
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

    def align_partial(self, primer: str, super_sequence: str, search_interval: tuple[int, int]) -> MatchResultDTO:
        """
        Align the primer sequence to the super sequence. You can set a search interval within the super sequence.
        :param primer: A primer sequence.
        :param super_sequence: A (longer) DNA sequence.
        :param search_interval: The search interval within the super sequence.
        :return: A MatchResultDTO with the alignment result or an empty match if an error occurs.
        """
        try:
            # Validate inputs
            if primer is None or super_sequence is None:
                logger.error(f"Invalid input to align_partial: primer={primer}, super_sequence length={len(super_sequence) if super_sequence is not None else None}")
                return MatchResultDTO(0, "", 0, 0, primer if primer is not None else "")
                
            if not isinstance(search_interval, tuple) or len(search_interval) != 2:
                logger.error(f"Invalid search interval format: {search_interval}")
                return MatchResultDTO(0, "", 0, 0, primer)
                
            # Ensure search interval is valid
            start, end = search_interval
            if start < 0:
                logger.warning(f"Negative start index in search interval: {start}, setting to 0")
                start = 0
                
            #if end > len(super_sequence):
            #    logger.warning(f"End index {end} exceeds sequence length {len(super_sequence)}, adjusting")
            #    end = len(super_sequence)
                
            if start >= end:
                logger.error(f"Invalid search interval: start {start} >= end {end}")
                return MatchResultDTO(0, "", 0, 0, primer)
            
            # Perform alignment
            match = self.align(
                primer_sequence=primer,
                dna_sequence=super_sequence[start:end],
                are_ends_eligible_for_bonus=(start == 0, end == len(super_sequence))
            )
            
            # Adjust indices to account for the search interval
            match.start_index += start
            match.end_index += start
            return match
            
        except Exception as e:
            logger.error(f"Error in align_partial: {str(e)}")
            return MatchResultDTO(0, "", 0, 0, primer if primer is not None else "")


    def align(
            self,
            primer_sequence: str,
            dna_sequence: str,
            are_ends_eligible_for_bonus: tuple[bool, bool]
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
        try:
            # Validate inputs
            if primer_sequence is None:
                logger.error(f"Invalid input to align: primer_sequence={primer_sequence}, dna_sequence length={len(dna_sequence) if dna_sequence is not None else None}")
                return MatchResultDTO(0, "", 0, 0, primer_sequence if primer_sequence is not None else "")

            if not primer_sequence:
                logger.error(f"Empty sequence provided: primer_sequence length={len(primer_sequence)}, dna_sequence length={len(dna_sequence)}")
                return MatchResultDTO(0, "", 0, 0, primer_sequence)

            if dna_sequence is None or not dna_sequence:
                # fail silent, this happens.
                return MatchResultDTO(0, "", 0, 0, primer_sequence)

            if not isinstance(are_ends_eligible_for_bonus, tuple) or len(are_ends_eligible_for_bonus) != 2:
                logger.warning(f"Invalid are_ends_eligible_for_bonus format: {are_ends_eligible_for_bonus}, using (False, False)")
                are_ends_eligible_for_bonus = (False, False)
            
            # construction
            try:
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
                        try:
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
                        except IndexError as e:
                            logger.error(f"Index error in alignment matrix calculation: {str(e)}, i={i}, j={j}")
                            continue
                        except Exception as e:
                            logger.error(f"Error in alignment calculation: {str(e)}")
                            continue

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
                    try:
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
                        else:
                            logger.warning(f"Unknown direction in traceback: {direction}")
                            break
                    except IndexError as e:
                        logger.error(f"Index error during traceback: {str(e)}, i={i}, j={j}")
                        break
                    except Exception as e:
                        logger.error(f"Error during traceback: {str(e)}")
                        break

                aligned_read = ''.join(reversed(aligned_read))

                # There were 5 in 2.3M sequences, which managed to set their index to -1 instead of 0.
                # Since their alignment worked otherwise, I set the index to min 0. Maybe their score is also off by 1 though.
                # The same for some alignments that managed to end outside the sequence.
                alignment_start_index_in_read = max(0, j - 2)
                alignment_end_index = min(alignment_start_index_in_read + len(aligned_read), len(dna_sequence))


                return MatchResultDTO(max_score, aligned_read, alignment_start_index_in_read, alignment_end_index, primer_sequence)
            except Exception as e:
                logger.error(f"Error in Smith-Waterman algorithm: {str(e)}")
                return MatchResultDTO(0, "", 0, 0, primer_sequence)
                
        except Exception as e:
            logger.error(f"Unexpected error in align method: {str(e)}")
            return MatchResultDTO(0, "", 0, 0, primer_sequence if primer_sequence is not None else "")
