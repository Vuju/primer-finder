from dataclasses import dataclass

@dataclass
class MatchResultDTO:
    """
    A dataclass representing a match result.
    """
    score: float = 0.0
    read: str = ''
    start_index: int = -1
    end_index: int = -1
    primer_sequence: str = ''
    quality_cutoff: float = 0.7

    def is_mismatch(self):
        """
        Checks if the match result is mismatch.
        :return: Boolean whether the match result is mismatch.
        """
        return True if self.start_index == -1 else False

