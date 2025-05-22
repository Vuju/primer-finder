from abc import ABC, abstractmethod
from typing import Generator, Tuple, Any

from primer_finder.matching.dtos.match_result_dto import MatchResultDTO


class Connector(ABC):
    """
    An abstract base class that defines the interface for all connectors for the primer finder class.
    """
    @abstractmethod
    def read_sequences(self, forward_primer, backward_primer, batch_size) -> (
            Generator)[tuple[Any, Any, MatchResultDTO, MatchResultDTO], Any, None]:
        """
        Returns a generator that yields each sequence in the input.
        """
        pass

    @abstractmethod
    def write_output(self, lock, information):
        """
        writes all the important information back to the output.
        information should be a list of a big tuple like:
        [(read_metadata, forward_match, backward_match, inter_primer_sequence, possible_orf)]
        """
        pass

    @abstractmethod
    def get_number_of_sequences(self) -> int:
        """
        Returns the number of sequences in the input.
        """
        pass