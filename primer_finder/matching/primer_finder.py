import logging
import time
from functools import partial
from multiprocessing import Pool
from typing import Generator, Any

from tqdm import tqdm

from primer_finder.config.config_loader import get_config_loader, get_search_parameter_objects
from primer_finder.orf.finder import list_possible_orf
from primer_finder.connectors.base import Connector
from primer_finder.matching.dtos.match_result_dto import MatchResultDTO
from primer_finder.matching.dtos.search_parameter_object import SearchParameterObject
from primer_finder.matching.regex import find_exact_match
from primer_finder.matching.smith_waterman import SmithWaterman

logger = logging.getLogger(__name__)

def chunker(iterable, sub_chunk_size: int) -> Generator[list, Any, None]:
    """
    Helper method that converts an itarable into a generator of chunks.
    :param iterable: list to be chunked.
    :param sub_chunk_size: size of each chunk.
    :return: Generator of chunks.
    """
    chunk = []
    for item in iterable:
        chunk.append(item)
        if len(chunk) >= sub_chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk

class PrimerFinder:
    """
    A configured instance of PrimerFinder.
    """
    def __init__(self,
                 connector: Connector,
                 custom_num_threads: int = None,
                 chunk_size: int = None,
                 search_area: float = None,
                 search_parameter_groups: list[SearchParameterObject] = None,
                 ):
        """
        A configured instance of PrimerFinder.
        :param connector: a connector object to handle I/O operations.
        :param custom_num_threads: Number of threads to use for parallel processing.
        :param chunk_size: Size of each chunk to be processed in parallel.
        :param search_area: Percentage of the distance between two primers used to reduce the search area for the second primer.
        :param search_parameter_groups: List of SearchParameterObjects used as query parameters.
        """
        config = get_config_loader().get_config()

        self.search_area = search_area or config["algorithm"]["search_area"]
        self.custom_num_threads = custom_num_threads or config["parallelization"]["num_threads"]
        self.chunk_size = chunk_size or config["parallelization"]["chunk_size"]
        self.buffer_flush_threshold = config["database"]["database_batch_size"] # todo: this should not be here
        self.search_parameter_groups = search_parameter_groups or get_search_parameter_objects()
        self.connector = connector
        self.smith_waterman = SmithWaterman()

    def find_all_primers(self):
        """
        Executes the primer finder algorithm with the set configuration.
        """
        try:
            logger.info("Getting the number of sequences.")
            _total_number_of_sequences = self.connector.get_number_of_sequences()
            for i, primer_datum in enumerate(self.search_parameter_groups):
                try:
                    sequences = self.connector.read_sequences(primer_datum.forward_primer, primer_datum.reverse_primer)
                    logger.info(f"Searching input sequences for primer pair {i + 1}.")
                    pbar = tqdm(total=_total_number_of_sequences)
                    worker = partial(self._process_sequences_chunk, primer_datum)

                    sequence_chunks = chunker(sequences, self.chunk_size)
                    with Pool(processes=self.custom_num_threads) as pool:
                        results_buffer = []
                        try:
                            for wbv in pool.imap(worker, sequence_chunks, chunksize=self.chunk_size//10):
                                results_buffer.extend(wbv)
                                # todo the buffering should be the responsibility of a connector, even though it makes connectors harder to write
                                if len(results_buffer) >= self.buffer_flush_threshold:
                                    try:
                                        if self.connector.write_output(results_buffer):
                                            results_buffer = []
                                    except Exception as e:
                                        logger.error(f"Error writing output batch: {str(e)}")
                                pbar.update(self.chunk_size)
                            if results_buffer:
                                retry_count = 0
                                max_retries = 5
                                while not self.connector.write_output(results_buffer) and retry_count < max_retries:
                                    time.sleep(5)
                                    retry_count += 1
                                if retry_count >= max_retries:
                                    logger.error("Failed to write remaining results after maximum retries")
                        except Exception as e:
                            logger.error(f"Error during sequence processing: {str(e)}")
                        finally:
                            pbar.close()
                except Exception as e:
                    logger.error(f"Error processing primer pair {i + 1}: {str(e)}")
        except Exception as e:
            logger.error(f"Fatal error in primer finding process: {str(e)}")
            raise

    def _process_sequences_chunk(self, query: SearchParameterObject, sequence_list):
        writeback_values = []
        for sequence_object in sequence_list:
            try:
                writeback_values.append(self._process_sequence(query, sequence_object))
            except Exception as e:
                logger.error(f"Error processing sequence {sequence_object[0] if len(sequence_object) > 0 else 'unknown'}: {str(e)}")
                # Add a placeholder result to maintain sequence count
                if len(sequence_object) > 0:
                    writeback_values.append((sequence_object[0], MatchResultDTO(), MatchResultDTO(), None, [], query.distance))
        return writeback_values

    def _process_sequence(self, query: SearchParameterObject, sequence_object):
        try:
            _sequence_found = False
            _orf_calculated = 0

            offset = int(query.distance * self.search_area)
            distance = query.distance
            def __limit_to_interval_after(i: int):
                return (i + distance - offset,
                        i + distance + len(query.reverse_primer) + offset)
            def __limit_to_interval_before(i: int):
                return (max(0, i - distance - len(query.forward_primer) - offset),
                        max(0, i - distance + offset))

            try:
                sequence_metadata, dna_sequence, forward_match, reverse_match = sequence_object
            except ValueError as e:
                logger.error(f"Invalid sequence object format: {str(e)}")
                return sequence_object[0] if len(sequence_object) > 0 else None, MatchResultDTO(), MatchResultDTO(), None, [], query.distance

            if dna_sequence is None:
                return sequence_metadata, MatchResultDTO(), MatchResultDTO(), None, [], query.distance

            dna_sequence = dna_sequence.strip()
            forward_search_interval, reverse_search_interval = (0, len(dna_sequence)), (0, len(dna_sequence))

            ## first check for exact matches
            try:
                if forward_match.is_mismatch():
                    forward_match = self._compute_regex_match(query.forward_primer,
                                                            query.forward_primer_regex,
                                                            dna_sequence)
            except Exception as e:
                logger.error(f"Error in forward primer regex matching: {str(e)}")
                forward_match = MatchResultDTO(0, "", 0, 0, query.forward_primer)

            try:
                if reverse_match.is_mismatch():
                    if not forward_match.is_mismatch():
                        reverse_search_interval = __limit_to_interval_after(forward_match.end_index)

                    reverse_match = self._compute_regex_match(query.reverse_primer,
                                                           query.reverse_primer_regex,
                                                           dna_sequence[reverse_search_interval[0]:reverse_search_interval[1]])
                    if not reverse_match.is_mismatch():
                        reverse_match.start_index += reverse_search_interval[0]
                        reverse_match.end_index += reverse_search_interval[0]
                        forward_search_interval = __limit_to_interval_before(reverse_match.start_index)
            except Exception as e:
                logger.error(f"Error in reverse primer regex matching: {str(e)}")
                reverse_match = MatchResultDTO(0, "", 0, 0, query.reverse_primer)

            ## for each missing exact match, try smith waterman:
            try:
                if forward_match.is_mismatch():
                    forward_match = self.smith_waterman.align_partial(
                        primer=query.forward_primer,
                        super_sequence=dna_sequence,
                        search_interval=forward_search_interval
                    )
                    score_threshold = (len(query.forward_primer)
                                    * self.smith_waterman.match_value
                                    * query.forward_cutoff)

                    if reverse_match.is_mismatch() and (forward_match.score > score_threshold):
                        reverse_search_interval = __limit_to_interval_after(forward_match.end_index)
            except Exception as e:
                logger.error(f"Error in forward primer Smith-Waterman alignment: {str(e)}")
                forward_match = MatchResultDTO(0, "", 0, 0, query.forward_primer)

            try:
                if reverse_match.is_mismatch():
                    reverse_match = self.smith_waterman.align_partial(
                        primer=query.reverse_primer,
                        super_sequence=dna_sequence,
                        search_interval=reverse_search_interval
                    )
            except Exception as e:
                logger.error(f"Error in reverse primer Smith-Waterman alignment: {str(e)}")
                reverse_match = MatchResultDTO(0, "", 0, 0, query.reverse_primer)

            ## work on getting the orf
            try:
                # Ensure indices are valid
                if forward_match.end_index < 0:
                    forward_match.end_index = 0
                if reverse_match.start_index < 0:
                    reverse_match.start_index = 0
                if forward_match.end_index > len(dna_sequence):
                    forward_match.end_index = len(dna_sequence)
                if reverse_match.start_index > len(dna_sequence):
                    reverse_match.start_index = len(dna_sequence)
                
                # Ensure end_index is not greater than start_index
                if forward_match.end_index > reverse_match.start_index:
                    inter_primer_region = ""
                else:
                    inter_primer_region = dna_sequence[forward_match.end_index:reverse_match.start_index]
                
                _sequence_found = len(inter_primer_region.strip()) > 0

                if not _sequence_found:
                    if (forward_match.score / max(1, len(forward_match.primer_sequence))) < (reverse_match.score / max(1, len(reverse_match.primer_sequence))):
                        forward_match = MatchResultDTO(0, "", 0, 0, query.forward_primer)
                    else:
                        reverse_match = MatchResultDTO(0, "", 0, 0, query.reverse_primer)
                possible_orfs = []
                if _sequence_found:
                    try:
                        possible_orfs = list_possible_orf(inter_primer_region, translation_table=query.protein_translation_table)
                        possible_orfs = ([]) if len(possible_orfs) == 0 else possible_orfs
                    except Exception as e:
                        logger.error(f"Error finding possible ORFs: {str(e)}")
                        possible_orfs = []
                else:
                    inter_primer_region = None
            except Exception as e:
                logger.error(f"Error processing inter-primer region: {str(e)}")
                inter_primer_region = None
                possible_orfs = []

            forward_match.quality_cutoff = query.forward_cutoff
            reverse_match.quality_cutoff = query.reverse_cutoff

            return sequence_metadata, forward_match, reverse_match, inter_primer_region, possible_orfs, query.distance
        except Exception as e:
            logger.error(f"Unexpected error in _process_sequence: {str(e)}")
            # Return a safe default value
            if isinstance(sequence_object, tuple) and len(sequence_object) > 0:
                return sequence_object[0], MatchResultDTO(), MatchResultDTO(), None, [], query.distance
            else:
                return None, MatchResultDTO(), MatchResultDTO(), None, [], query.distance

    def _compute_regex_match(self, primer, primer_regex, read):
        try:
            score = 0
            read_match = ''
            
            if primer_regex is None or read is None:
                logger.warning(f"Invalid input for regex match: primer_regex={primer_regex}, read length={len(read) if read is not None else None}")
                return MatchResultDTO(0, "", -1, -1, primer)
                
            try:
                index, end_index = find_exact_match(primer_regex, read)
            except Exception as e:
                logger.error(f"Error in regex matching: {str(e)}")
                return MatchResultDTO(0, "", -1, -1, primer)
                
            if index != -1:
                score = len(primer) * self.smith_waterman.match_value
                try:
                    read_match = read[index:end_index]
                except IndexError as e:
                    logger.error(f"Index error when extracting match: {str(e)}, indices: {index}:{end_index}, read length: {len(read)}")
                    index, end_index = max(0, index), min(end_index, len(read))
                    read_match = read[index:end_index] if index < end_index else ""
                    
            return MatchResultDTO(score, read_match, index, end_index, primer)
        except Exception as e:
            logger.error(f"Unexpected error in _compute_regex_match: {str(e)}")
            return MatchResultDTO(0, "", -1, -1, primer)
