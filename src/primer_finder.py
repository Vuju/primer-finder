import argparse
import gzip
import logging

from functools import partial
from multiprocessing import Pool, Lock

import pandas as pd
from tqdm import tqdm
from typing import TextIO, Any

from src.match_result import MatchResult
from src.orf_finder import list_possible_orf, OrfFinder
from src.primer_data_dto import PrimerDataDTO, get_primer_dto_from_args
from src.primer_finder_regex import *
from src.smith_waterman import SmithWaterman

logger = logging.getLogger(__name__)
smith_waterman: SmithWaterman

# todo: create options.py creation + parameterization

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process sequence alignment parameters.")

    parser.add_argument("--primer_finder", type=bool, default=True,
                        help="Flag as false to disable the primer-searching algorithm.")
    parser.add_argument("--orf_matching", type=bool, default=False,
                        help="Flag as false to disable the orf-decision algorithm.")

    parser.add_argument("--search_area", type=float, default=0.2, help="This value will determine, "
                                                        "how much extra area the smith waterman algorithm will search, "
                                                        "if the other primer has already been found with enough certainty (set by '--sw_cutoff').")
    parser.add_argument("--smith_waterman_score_cutoff", type=float, default=0.8, help="Smith-Waterman score cutoff (default: 0.8)")

    parser.add_argument("--primer_information", type=str, default="./data/primer-information.csv",
                        help="CSV list of forward and reverse primer sequence, as well as the expected distance inbetween.")

    parser.add_argument("--muscle_path", type=str, default="/mnt/c/Users/Julian/bin/muscle",
                        help="Path to the muscle binary/executable. I run with version 5.3, and it will be used to run 'muscle_path -align tmp_in.fasta -out tmp_out.fasta'")

    parser.add_argument("--input_file_path", type=str, default="./data/DB.COX1.fna", help="Path to input sequence file")
    parser.add_argument("--output_file_path", type=str, default="./data/primer-finder-result.csv",
                        help="Path to output results file")
    parser.add_argument("--orf_matching_threshold", type=int, default=10,
                        help="Minimum number of similar sequences required to match an orf")
    parser.add_argument("--orf_matching_upper_threshold", type=int, default=50,
                        help="Limit of similar sequences used to match an orf")
    parser.add_argument("--protein_translation_table", type=Any, default=5,
                        help="Translation table for Bio.Seq translate(). This is used in orf_finder.")
    parser.add_argument("--num_threads", type=int, default=None,
                        help="Translation table for Bio.Seq translate(). This is used in orf_finder.")

    return parser.parse_args()


def compute_arguments():
    args = parse_arguments()

    ### hardcoded parameters:
    args.gap_penalty = -2
    args.gap3_penalty = -2
    args.end_of_read_bonus = 1
    args.substitution_function = None
    args.chunksize = 100
    args.e_value = 1000
    ### end hardcoded parameters

    args.primer_data = []

    with open(args.primer_information, "r") as primer_info_file:
        for line in primer_info_file.readlines():
            line_data = line.split(",")
            line_data = [s.strip() for s in line_data]

            f_primer_regex = regex_builder(line_data[0])
            b_primer_regex = regex_builder(line_data[1])

            entry = {
                "f_primer": line_data[0],
                "b_primer": line_data[1],
                "distance": line_data[2],
                "f_primer_regex": f_primer_regex,
                "b_primer_regex": b_primer_regex
            }
            args.primer_data.append(entry)


    return args



### Function definitions
def read_pairs(file_path):
    file: TextIO
    if file_path.endswith('.gz'):
        file = gzip.open(file_path, 'rt')
    else:
        file = open(file_path, 'r', encoding="UTF-8")

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


def init(l):
    global lock
    lock = l


def write_output_to_file(output_path, read_metadata, forward_match, backward_match, sequence, possible_orf):
    with lock, open(output_path, 'a') as out_file:
        out_file.write(read_metadata.replace('|', ';').replace(',', ';').strip() +
                       f"{forward_match.score};{forward_match.read};{forward_match.start_index};" +
                       f"{backward_match.score};{backward_match.read};{backward_match.start_index};{sequence};{str(possible_orf)}\n")


def compute_regex_match(primer, primer_regex, read):
    score = 0
    read_match = ''
    index, end_index = find_exact_match(primer_regex, read)
    if index != -1:
        score = len(primer) * smith_waterman.match_value
        read_match = read[index:end_index]
    return MatchResult(score, read_match, index, end_index)


### the main processing function

def process_sequence(primer_data: PrimerDataDTO, sequence_object):
    # initializing all flags
    _sequence_found = False
    _orf_calculated = 0
        
    
    offset = int(primer_data.distance * primer_data.search_area)
    distance = primer_data.distance

    def limit_to_interval_after(i: int):
        return (i + distance - offset,
                i + distance + len(primer_data.backward_primer) + offset)
    def limit_to_interval_before(i: int):
        return (max(0, i - distance - len(primer_data.forward_primer) - offset),
                max(0, i - distance + offset))

    sequence_metadata, dna_sequence = sequence_object
    dna_sequence = dna_sequence.strip()
    forward_search_interval, backward_search_interval = (0, len(dna_sequence)), (0, len(dna_sequence))

    ## first check for exact matches
    forward_match = compute_regex_match(primer_data.forward_primer,
                                  primer_data.forward_primer_regex,
                                  dna_sequence)
    if forward_match.start_index != -1:
        backward_search_interval = limit_to_interval_after(forward_match.end_index)

    backward_match = compute_regex_match(primer_data.backward_primer,
                                  primer_data.backward_primer_regex,
                                  dna_sequence[backward_search_interval[0]:backward_search_interval[1]])
    if backward_match.start_index != -1:
        backward_match.start_index += backward_search_interval[0]
        backward_match.end_index += backward_search_interval[0]
        forward_search_interval = limit_to_interval_before(backward_match.start_index)

    ## for each missing exact match, try smith waterman:
    if forward_match.is_mismatch():
        forward_match = smith_waterman.align_partial(
            primer=primer_data.forward_primer,
            super_sequence=dna_sequence,
            search_interval=forward_search_interval
        )

        score_threshold = (len(primer_data.forward_primer)
                           * smith_waterman.match_value
                           * primer_data.smith_waterman_score_cutoff)

        if backward_match.is_mismatch() and (forward_match.score > score_threshold):
            backward_search_interval = limit_to_interval_after(forward_match.end_index)

    if backward_match.is_mismatch():
        backward_match = smith_waterman.align_partial(
            primer=primer_data.backward_primer,
            super_sequence=dna_sequence,
            search_interval=backward_search_interval
        )

    ## work on getting the orf
    inter_primer_region = dna_sequence[forward_match.end_index:backward_match.start_index]
    _sequence_found = len(inter_primer_region.strip()) > 0
    
    if _sequence_found:
         
        possible_orf = list_possible_orf(inter_primer_region, translation_table=primer_data.translation_table)
        if len(possible_orf) == 0:
            _orf_calculated = -1
        elif len(possible_orf) == 1:
            _orf_calculated = 1
        elif len(possible_orf) >= 2:
            _orf_calculated = 2
        else:
            _orf_calculated = -2
        
        possible_orf = ([]) if _orf_calculated <= 0 else possible_orf

        write_output_to_file(primer_data.output_file_path, sequence_metadata, forward_match, backward_match, dna_sequence, possible_orf)


### main script

currentRead = ""


def get_number_of_sequences_in(input_file_path):
    count = 0
    file: TextIO
    if input_file_path.endswith('.gz'):
        file = gzip.open(input_file_path, 'rt')
    else:
        file = open(input_file_path, 'r', encoding="UTF-8")
    for line in file.readlines():
        if line.startswith('>'):
            count += 1
    file.close()
    return count


if __name__ == "__main__":
    logging.basicConfig(filename='primer-finder_log.log', level=logging.INFO)
    logging.getLogger().addHandler(logging.StreamHandler())

    args = compute_arguments()
    smith_waterman = SmithWaterman(
        gap_penalty=args.gap_penalty,
        triple_gap_penalty=args.gap3_penalty,
        end_of_read_bonus_value = args.end_of_read_bonus,
        sub_function=args.substitution_function
    )

    total_number_of_sequences = 0
    for i, primer_pair in enumerate(args.primer_data):

        primer_data = get_primer_dto_from_args(args, i)

        if i == 0:
            logger.info(f"Getting length of input file.")
            total_number_of_sequences=get_number_of_sequences_in(primer_data.input_file_path)

        if args.primer_finder:
            with open(primer_data.output_file_path, 'w') as output_file:
                output_file.write(
                    "BOLD ID;Read ID;Country;Phylum;Class;Order;Family;Genus;Species;f_score;f_match;f_index;b_score;b_match;b_index;read;possible_orfs\n"
                )

            pairs = read_pairs(primer_data.input_file_path)
            lock = Lock()
            logger.info(f"Searching input sequences for primer pair {i + 1}.")
            pbar = tqdm(total=total_number_of_sequences)
            worker = partial(process_sequence, primer_data)
            with Pool(processes=args.num_threads, initializer=init, initargs=(lock,)) as pool:
                for _ in pool.imap(worker, pairs, chunksize=args.chunksize):
                    pbar.update(1)
            pbar.close()

            logger.info(f"Primer Finder output has been written to {primer_data.output_file_path}")

        if args.orf_matching:
            logger.info(f"Starting orf-matching process.")
            orf_finder = OrfFinder(
                lower_reference_threshold=args.orf_matching_threshold,
                upper_reference_threshold=args.orf_matching_upper_threshold,
                translation_table=args.protein_translation_table,
                muscle_path=args.muscle_path,
                e_value_threshold=args.e_value
            )
            all_entries = pd.read_csv(primer_data.output_file_path, sep=";")
            solved = orf_finder.solve_orfs_for_df(df=all_entries)
            logger.info("Starting write-back.")
            solved.to_csv(primer_data.output_file_path)

            logger.info(f"Orf matching output has been written to {primer_data.output_file_path}")
