import gzip
import logging

from functools import partial
from multiprocessing import Pool, Lock

from tqdm import tqdm
from typing import TextIO

from src.config.constants import INPUT_FILE_PATH, OUTPUT_FILE_PATH, \
    PRIMER_INFORMATION_PATH, CUSTOM_NUM_THREADS, CHUNKSIZE, SEARCH_AREA, SMITH_WATERMAN_SCORE_CUTOFF, \
    PROTEIN_TRANSLATION_TABLE
from src.match_result import MatchResult
from src.orf_finder import list_possible_orf
from src.primer_data_dto import PrimerDataDTO, primer_info_from_string
from src.primer_finder_regex import *
from src.smith_waterman import SmithWaterman

logger = logging.getLogger(__name__)

lock = Lock()
def _init_lock(l):
    global lock
    lock = l

class PrimerFinder:

    def __init__(self,
                 smith_waterman: SmithWaterman = None,
                 input_file: str = INPUT_FILE_PATH,
                 output_file: str = OUTPUT_FILE_PATH,
                 primer_information_file: str = PRIMER_INFORMATION_PATH,
                 custom_num_threads: int = CUSTOM_NUM_THREADS,
                 chunk_size: int = CHUNKSIZE,
                 search_area: float = SEARCH_AREA,
                 smith_waterman_score_cutoff: float = SMITH_WATERMAN_SCORE_CUTOFF,
                 translation_table = PROTEIN_TRANSLATION_TABLE,
    ):
        self.primer_data = []
        self.smith_waterman = smith_waterman or SmithWaterman()
        self.input_file = input_file
        self.output_file = output_file
        self.primer_information_file = primer_information_file
        self.custom_num_threads = custom_num_threads
        self.chunk_size = chunk_size
        self.search_area = search_area
        self.smith_waterman_score_cutoff = smith_waterman_score_cutoff
        self.translation_table = translation_table

    def find_all_primers(self):

        self.get_primer_information()
        sequences = self.read_sequences()
        _lock = Lock()

        logger.info("Getting the number of sequences.")
        _total_number_of_sequences = self.count_number_of_sequences()
        with open(self.output_file, 'w') as output_file:
            output_file.write(
                "BOLD ID;Read ID;Country;Phylum;Class;Order;Family;Genus;Species;f_score;f_match;f_index;b_score;b_match;b_index;read;possible_orfs\n"
            )
        for i, primer_datum in enumerate(self.primer_data):
            logger.info(f"Searching input sequences for primer pair {i + 1}.")
            pbar = tqdm(total=_total_number_of_sequences)
            worker = partial(self.process_sequence, primer_datum)
            with Pool(processes=self.custom_num_threads, initializer=_init_lock, initargs=(_lock,)) as pool:
                for _ in pool.imap(worker, sequences, chunksize=self.chunk_size):
                    pbar.update(1)
            pbar.close()

    def get_primer_information(self):
        with open(self.primer_information_file, "r") as primer_info_file:
            for line in primer_info_file.readlines():
                entry = primer_info_from_string(line)
                self.primer_data.append(entry)



    def read_sequences(self):
        file: TextIO
        if self.input_file.endswith('.gz'):
            file = gzip.open(self.input_file, 'rt')
        else:
            file = open(self.input_file, 'r', encoding="UTF-8")

        while True:
            metadata_line = file.readline()
            sequence_lines = file.readline()
            if not metadata_line or not sequence_lines:
                break
            while True:
                pos = file.tell()
                line_next = file.readline()
                if not line_next:
                    break
                if line_next.strip() == "":
                    break
                if line_next.startswith('>'):
                    file.seek(pos)
                    break
                else:
                    sequence_lines += line_next.strip()

            yield metadata_line, sequence_lines

    def compute_regex_match(self, primer, primer_regex, read):
        score = 0
        read_match = ''
        index, end_index = find_exact_match(primer_regex, read)
        if index != -1:
            score = len(primer) * self.smith_waterman.match_value
            read_match = read[index:end_index]
        return MatchResult(score, read_match, index, end_index)

    def process_sequence(self, query: PrimerDataDTO, sequence_object):
        _sequence_found = False
        _orf_calculated = 0

        offset = int(query.distance * self.search_area)
        distance = query.distance
        def limit_to_interval_after(i: int):
            return (i + distance - offset,
                    i + distance + len(query.backward_primer) + offset)
        def limit_to_interval_before(i: int):
            return (max(0, i - distance - len(query.forward_primer) - offset),
                    max(0, i - distance + offset))

        sequence_metadata, dna_sequence = sequence_object
        dna_sequence = dna_sequence.strip()
        forward_search_interval, backward_search_interval = (0, len(dna_sequence)), (0, len(dna_sequence))

        ## first check for exact matches
        forward_match = self.compute_regex_match(query.forward_primer,
                                                 query.forward_primer_regex,
                                                 dna_sequence)
        if forward_match.start_index != -1:
            backward_search_interval = limit_to_interval_after(forward_match.end_index)

        backward_match = self.compute_regex_match(query.backward_primer,
                                                  query.backward_primer_regex,
                                      dna_sequence[backward_search_interval[0]:backward_search_interval[1]])
        if backward_match.start_index != -1:
            backward_match.start_index += backward_search_interval[0]
            backward_match.end_index += backward_search_interval[0]
            forward_search_interval = limit_to_interval_before(backward_match.start_index)

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
                backward_search_interval = limit_to_interval_after(forward_match.end_index)

        if backward_match.is_mismatch():
            backward_match = self.smith_waterman.align_partial(
                primer=query.backward_primer,
                super_sequence=dna_sequence,
                search_interval=backward_search_interval
            )

        ## work on getting the orf
        inter_primer_region = dna_sequence[forward_match.end_index:backward_match.start_index]
        _sequence_found = len(inter_primer_region.strip()) > 0

        if _sequence_found:
            possible_orfs = list_possible_orf(inter_primer_region, translation_table=self.translation_table)
            possible_orfs = ([]) if len(possible_orfs) == 0 else possible_orfs

            self.write_output_to_file(sequence_metadata, forward_match, backward_match, inter_primer_region, possible_orfs)


    def write_output_to_file(self, read_metadata, forward_match, backward_match, inter_primer_sequence, possible_orf):
        with lock, open(self.output_file, 'a') as out_file:
            out_file.write(read_metadata.replace('|', ';').replace(',', ';').strip()
                           + f"{forward_match.score};{forward_match.read};{forward_match.start_index};"
                           + f"{backward_match.score};{backward_match.read};{backward_match.start_index};"
                           + f"{inter_primer_sequence};{str(possible_orf)}\n")

    def count_number_of_sequences(self):
        count = 0
        file: TextIO
        if self.input_file.endswith('.gz'):
            file = gzip.open(self.input_file, 'rt')
        else:
            file = open(self.input_file, 'r', encoding="UTF-8")
        for line in file.readlines():
            if line.startswith('>'):
                count += 1
        file.close()
        return count
