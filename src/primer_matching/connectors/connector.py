from abc import ABC, abstractmethod
from typing import Generator, Tuple, Any


class Connector(ABC):
    """
    An abstract base class that defines the interface for all connectors for the primer finder class.
    """
    @abstractmethod
    def read_sequences(self) -> Generator[Tuple[str, str], Any, None]:
        """
        Returns a generator that yields each sequence in the input.
        """
        pass

    @abstractmethod
    def write_output(self, read_metadata, forward_match, backward_match, inter_primer_sequence, possible_orf):
        """
        writes all the important information back to the output.
        """
        pass

    @abstractmethod
    def get_number_of_sequences(self) -> int:
        """
        Returns the number of sequences in the input.
        """
        pass