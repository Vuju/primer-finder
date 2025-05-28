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


    #------------------------------------ ORF matching methods ---------------------------------
    @abstractmethod
    def read_pairs_chunk(self, chunk_size):
        """
        Reads a chunk of pairs of sequences from the primer-pairs table.
        :param chunk_size:
        :return:
        """
        pass

    @abstractmethod
    def write_pair_chunk(self, solved):
        """
        Updates the primer-pairs table with the solved pairs.
        :param solved:
        :return:
        """
        pass

    @abstractmethod
    def get_remaining_unsolved_count(self):
        """
        Returns the number of unsolved entries in the primer-pairs table.
        :return:
        """
        pass

    @abstractmethod
    def get_next_unsolved_sequence(self):
        """
        Fetches a single unsolved entry from the primer-pairs table.
        :return:
        """
        pass

    @abstractmethod
    def fetch_sampled_solved_related_sequences(self, current_entry, level, lower_reference_threshold, upper_reference_threshold):
        """
        Returns a comparison group for the given entry, based on taxonomy.
        :param current_entry:
        :param level:
        :param lower_reference_threshold:
        :param upper_reference_threshold:
        :return:
        """
        pass

    @abstractmethod
    def fetch_unsolved_related_sequences(self, current_entry, level):
        """
        Fetches related sequences for the given entry that are unsolved as well.
        :param current_entry:
        :param level:
        :return:
        """
        pass
