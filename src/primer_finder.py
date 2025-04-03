import argparse
import gzip
import logging

from functools import partial
from multiprocessing import Pool, Lock

import pandas as pd
from tqdm import tqdm
from typing import TextIO, Any

from src.match_result import MatchResult
from src.orf_finder import list_possible_orf, solve_orfs_for_df
from src.primer_data_dto import PrimerDataDTO, get_primer_dto_from_args
from src.primer_finder_regex import *
from src.smith_waterman import smith_waterman

logger = logging.getLogger(__name__)

# optional: improve offsets to be more accurate
# optional: create options.ini creation + parameterization?




def parse_arguments():
    parser = argparse.ArgumentParser(description="Process sequence alignment parameters.")

    parser.add_argument("--primer_finder", type=bool, default=True,
                        help="Flag as false to disable the primer-searching algorithm.")
    parser.add_argument("--orf_matching", type=bool, default=True,
                        help="Flag as false to disable the orf-decision algorithm.")

    parser.add_argument("--search_area", type=float, default=0.2, help="This value will determine, "
                                                        "how much extra area the smith waterman algorithm will search, "
                                                        "if the other primer has already been found with enough certainty (set by '--sw_cutoff').")
    parser.add_argument("--sw_score_cutoff", type=float, default=0.8, help="Smith-Waterman score cutoff (default: 0.8)")

    parser.add_argument("--primer_information", type=str, default="./data/primer-information.csv",
                        help="CSV list of forward and reverse primer sequence, as well as the expected distance inbetween.")

    parser.add_argument("--input_file_path", type=str, default="./data/DB.COX1.fna", help="Path to input sequence file")
    parser.add_argument("--output_file_path", type=str, default="./data/primer-finder-result.csv",
                        help="Path to output results file")
    parser.add_argument("--orf_matching_threshold", type=int, default=4,
                        help="Minimum number of similar sequences required to match an orf")
    parser.add_argument("--orf_matching_upper_threshold", type=int, default=50,
                        help="Limit of similar sequences used to match an orf")
    parser.add_argument("--protein_translation_table", type=Any, default=5,
                        help="Translation table for Bio.Seq translate(). This is used in orf_finder.")

    return parser.parse_args()


def compute_arguments():
    args = parse_arguments()

    ### hardcoded parameters:
    # substitution function (see later)
    args.end_of_read_bonus = 1
    args.gap_penalty = -2
    args.gap3_penalty = -2
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

def substitution_function(p, r) -> float:
    match r:
        case 'A':
            return 2 if p in "AWMRDHVN" else -1
        case 'C':
            return 2 if p in "CSMYBHVN" else -1
        case 'G':
            return 2 if p in "GSKRBDVN" else -1
        case 'T':
            return 2 if p in "TWKYBDHN" else -1
        case default:
            raise Exception(f"unknown literal in read sequence: '{r}'")


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
        score = len(primer) * substitution_function('A', 'A')
        read_match = read[index:end_index]
    return MatchResult(score, read_match, index, end_index)


def compute_smith_waterman(primer, read, skip, skip3, substitution, end_bonus):
    score, _, read_match, index = smith_waterman(
        primer=primer,
        read=read,
        gap=skip,
        gap3=skip3,
        substitution_function=substitution,
        end_of_read_bonus=end_bonus
    )
    return MatchResult(score, read_match, index, index + len(read_match))


### the main processing function

def process_pair(primer_data: PrimerDataDTO, pair):
    # initializing all flags
    _sequence_found = False
    _orf_calculated = 0
        
    
    offset = int(primer_data.distance * primer_data.search_area)
    distance = primer_data.distance

    read_metadata, sequence = pair
    read = sequence.strip()
    f_search_interval, b_search_interval = (0, len(read)), (0, len(read))

    ## first check for exact matches
    f_match = compute_regex_match(primer_data.f_primer, primer_data.f_primer_regex, read)
    if f_match.start_index != -1:
        b_search_interval = (f_match.end_index + distance - offset, f_match.end_index + distance + len(primer_data.b_primer) + offset)

    b_match = compute_regex_match(primer_data.b_primer, primer_data.b_primer_regex, read[b_search_interval[0]:b_search_interval[1]])
    if b_match.start_index != -1:
        b_match.start_index += b_search_interval[0]
        b_match.end_index += b_search_interval[0]
        f_search_interval = (max(0, b_match.start_index - distance + offset), max(0, b_match.start_index - distance - len(primer_data.f_primer) + offset))


    ## for each missing exact match, try smith waterman:
    if f_match.start_index == -1:
        f_match = compute_smith_waterman(
            primer=primer_data.f_primer,
            read=read[f_search_interval[0]:f_search_interval[1]],
            skip=primer_data.sw_gap,
            skip3=primer_data.sw_gap3,
            substitution=substitution_function,
            end_bonus=args.end_of_read_bonus
        )
        f_match.start_index += f_search_interval[0]
        f_match.end_index += f_search_interval[0]

        score_threshold = len(primer_data.f_primer) * substitution_function('A', 'A') * primer_data.sw_score_cutoff

        if (b_match.start_index == -1) and (f_match.score > score_threshold):
            b_search_interval = (f_match.end_index + distance - offset, f_match.end_index + distance + len(primer_data.b_primer) + offset)

    if b_match.start_index == -1:
        b_match = compute_smith_waterman(
            primer=primer_data.b_primer,
            read=read[b_search_interval[0]:b_search_interval[1]],
            skip=primer_data.sw_gap,
            skip3=primer_data.sw_gap3,
            substitution=substitution_function,
            end_bonus=args.end_of_read_bonus
        )
        b_match.start_index += b_search_interval[0]
        b_match.end_index += b_search_interval[0]

    ## work on getting the orf
    sequence = read[f_match.end_index:b_match.start_index]
    _sequence_found = len(sequence.strip()) > 0
    
    if _sequence_found:
         
        possible_orf = list_possible_orf(sequence, translation_table=primer_data.translation_table)
        if len(possible_orf) == 0:
            _orf_calculated = -1
        elif len(possible_orf) == 1:
            _orf_calculated = 1
        elif len(possible_orf) >= 2:
            _orf_calculated = 2
        else:
            _orf_calculated = -2
        
        possible_orf = ([]) if _orf_calculated <= 0 else possible_orf

        write_output_to_file(primer_data.output_file_path, read_metadata, f_match, b_match, sequence, possible_orf)


### main script

currentRead = ""


def get_number_of_sequences_in(input_file_path):
    count = 0
    with open(input_file_path, 'r') as input_file:
        for line in input_file.readlines():
            if line.startswith('>'):
                count += 1
    return count


if __name__ == "__main__":
    logging.basicConfig(filename='primer-finder_log.log', level=logging.INFO)
    logging.getLogger().addHandler(logging.StreamHandler())

    args = compute_arguments()
    for i, primer_pair in enumerate(args.primer_data):
        primer_data = get_primer_dto_from_args(args, i)

        if args.primer_finder:
            with open(primer_data.output_file_path, 'w') as output_file:
                output_file.write(
                    "BOLD ID;Read ID;Country;Phylum;Class;Order;Family;Genus;Species;f_score;f_match;f_index;b_score;b_match;b_index;read;possible_orfs\n"
                )

            pairs = read_pairs(primer_data.input_file_path)
            lock = Lock()

            logger.info(f"Getting length of input {i+1}.")
            pbar = tqdm(total=get_number_of_sequences_in(primer_data.input_file_path))
            worker = partial(process_pair, primer_data)
            with Pool(initializer=init, initargs=(lock,)) as pool:
                for _ in pool.imap(worker, pairs, chunksize=args.chunksize):
                    pbar.update(1)
            pbar.close()

            logger.info(f"Primer Finder output has been written to {primer_data.output_file_path}")

        if args.orf_matching:
            logger.info(f"Starting orf-matching process.")
            df = pd.read_csv(primer_data.output_file_path, sep=";")
            solved = solve_orfs_for_df(
                df=df,
                threshold=args.orf_matching_threshold,
                upper_threshold=args.orf_matching_upper_threshold,
                translation_table=args.protein_translation_table,
                e_value=args.e_value
            )
            solved.to_csv(primer_data.output_file_path)

            logger.info(f"Orf matching output has been written to {primer_data.output_file_path}")
