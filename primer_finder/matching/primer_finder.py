import logging
import time
from functools import partial
from multiprocessing import Pool, Lock

from tqdm import tqdm

from primer_finder.config.config_loader import get_config_loader, get_primer_information
from primer_finder.orf.finder import list_possible_orf
from primer_finder.connectors.base import Connector
from primer_finder.matching.dtos.match_result_dto import MatchResultDTO
from primer_finder.matching.dtos.primer_data_dto import PrimerDataDTO
from primer_finder.matching.regex import find_exact_match
from primer_finder.matching.smith_waterman import SmithWaterman

logger = logging.getLogger(__name__)

lock = Lock()
def _init_lock(l):
    global lock
    lock = l

def chunker(iterable, sub_chunk_size):
    """Collect items from an iterable into chunks of size sub_chunksize"""
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
                 smith_waterman: SmithWaterman = None,
                 primer_information_file: str = None,
                 custom_num_threads: int = None,
                 chunk_size: int = None,
                 search_area: float = None,
                 smith_waterman_score_cutoff: float = None,
                 translation_table = None,
                 ):
        config = get_config_loader().get_config()

        self.primer_information_file = primer_information_file or config["paths"]["primer_information"]
        self.search_area = search_area or config["algorithm"]["search_area"]
        self.smith_waterman_score_cutoff = smith_waterman_score_cutoff or config["algorithm"]["smith_waterman_score_cutoff"]
        self.translation_table = translation_table or config["algorithm"]["protein_translation_table"]
        self.custom_num_threads = custom_num_threads or config["parallelization"]["num_threads"]
        self.chunk_size = chunk_size or config["parallelization"]["chunk_size"]
        self.db_batch_size = config["parallelization"]["database_batch_size"]
        self.primer_data = []
        self.connector = connector
        self.smith_waterman = smith_waterman or SmithWaterman()

    def find_all_primers(self):
        """
        Executes the primer finder algorithm with the set configuration.
        """
        get_primer_information(self.primer_information_file, self.primer_data)
        _lock = Lock()

        logger.info("Getting the number of sequences.")
        _total_number_of_sequences = self.connector.get_number_of_sequences()
        for i, primer_datum in enumerate(self.primer_data):
            sequences = self.connector.read_sequences(primer_datum.forward_primer, primer_datum.backward_primer)
            logger.info(f"Searching input sequences for primer pair {i + 1}.")
            pbar = tqdm(total=_total_number_of_sequences)
            worker = partial(self._process_sequences_chunk, primer_datum)

            sequence_chunks = chunker(sequences, self.chunk_size)
            with Pool(processes=self.custom_num_threads) as pool:
                results_buffer = []
                for wbv in pool.imap(worker, sequence_chunks, chunksize=self.chunk_size//10):
                    results_buffer.extend(wbv)
                    if len(results_buffer) >= self.db_batch_size:  # Adjust batch size
                        if self.connector.write_output("", results_buffer):
                            results_buffer = []
                    pbar.update(self.chunk_size)
                if results_buffer:
                    while not self.connector.write_output("", results_buffer):
                        time.sleep(5)
            pbar.close()

    def _process_sequences_chunk(self, query: PrimerDataDTO, sequence_list):
        writeback_values = []
        for sequence_object in sequence_list:
            writeback_values.append(self._process_sequence(query, sequence_object))
        return writeback_values

    def _process_sequence(self, query: PrimerDataDTO, sequence_object):
        _sequence_found = False
        _orf_calculated = 0

        offset = int(query.distance * self.search_area)
        distance = query.distance
        def __limit_to_interval_after(i: int):
            return (i + distance - offset,
                    i + distance + len(query.backward_primer) + offset)
        def __limit_to_interval_before(i: int):
            return (max(0, i - distance - len(query.forward_primer) - offset),
                    max(0, i - distance + offset))

        sequence_metadata, dna_sequence, forward_match, backward_match = sequence_object

        if dna_sequence is None:
            return sequence_metadata, MatchResultDTO(), MatchResultDTO(), None, []
        dna_sequence = dna_sequence.strip()
        forward_search_interval, backward_search_interval = (0, len(dna_sequence)), (0, len(dna_sequence))

        ## first check for exact matches
        if forward_match.is_mismatch():
            forward_match = self._compute_regex_match(query.forward_primer,
                                                      query.forward_primer_regex,
                                                      dna_sequence)
        if backward_match.is_mismatch():
            if not forward_match.is_mismatch():
                backward_search_interval = __limit_to_interval_after(forward_match.end_index)

            backward_match = self._compute_regex_match(query.backward_primer,
                                                       query.backward_primer_regex,
                                                       dna_sequence[backward_search_interval[0]:backward_search_interval[1]])
            if not backward_match.is_mismatch():
                backward_match.start_index += backward_search_interval[0]
                backward_match.end_index += backward_search_interval[0]
                forward_search_interval = __limit_to_interval_before(backward_match.start_index)

        ## for each missing exact match, try smith waterman:
        if forward_match.is_mismatch():
            forward_match = self.smith_waterman.align_partial(
                primer=query.forward_primer,
                super_sequence=dna_sequence,
                search_interval=forward_search_interval
            )
            score_threshold = (len(query.forward_primer)
                               * self.smith_waterman.match_value
                               * self.smith_waterman_score_cutoff)

            if backward_match.is_mismatch() and (forward_match.score > score_threshold):
                backward_search_interval = __limit_to_interval_after(forward_match.end_index)

        if backward_match.is_mismatch():
            backward_match = self.smith_waterman.align_partial(
                primer=query.backward_primer,
                super_sequence=dna_sequence,
                search_interval=backward_search_interval
            )

        ## work on getting the orf
        inter_primer_region = dna_sequence[forward_match.end_index:backward_match.start_index]
        _sequence_found = len(inter_primer_region.strip()) > 0

        possible_orfs = []
        if _sequence_found:
            possible_orfs = list_possible_orf(inter_primer_region, translation_table=self.translation_table)
            possible_orfs = ([]) if len(possible_orfs) == 0 else possible_orfs

        return sequence_metadata, forward_match, backward_match, inter_primer_region, possible_orfs

    def _compute_regex_match(self, primer, primer_regex, read):
        score = 0
        read_match = ''
        index, end_index = find_exact_match(primer_regex, read)
        if index != -1:
            score = len(primer) * self.smith_waterman.match_value
            read_match = read[index:end_index]
        return MatchResultDTO(score, read_match, index, end_index, primer)
