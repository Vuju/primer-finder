from abc import ABC, abstractmethod
from typing import Generator, Any

from dtos.match_result_dto import MatchResultDTO


class Connector(ABC):
    """
    An abstract base class that defines the interface for all connectors for the primer finder class.
    """
    @abstractmethod
    def read_sequences_for_matching(self, forward_primer, backward_primer, batch_size) -> (
            Generator)[tuple[Any, Any, MatchResultDTO, MatchResultDTO], Any, None]:
        """
        Returns a generator that yields each sequence in the input.
        """
        pass

    @abstractmethod
    def write_output_matches(self, information):
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


    #------------------------------------ ORF matching methods ---------------------------------
    @abstractmethod
    def read_pairs_chunk(self, chunk_size):
        pass

    @abstractmethod
    def write_pair_chunk(self, solved):
        pass

    @abstractmethod
    def get_remaining_unsolved_count(self):
        pass

    @abstractmethod
    def get_next_unsolved_sequence(self):
        pass

    @abstractmethod
    def find_comparison_group(self, current_entry, level, lower_reference_threshold, upper_reference_threshold):
        pass

    @abstractmethod
    def fetch_related_sequences(self, current_entry, level):
        pass

    @abstractmethod
    def update_sequences_with_results(self, solved):
        pass