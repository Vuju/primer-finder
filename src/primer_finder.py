import argparse
import gzip
import pickle
from functools import partial
from multiprocessing import Pool, Lock
from typing import TextIO

from src.smith_waterman import smith_waterman
from primer_finder_regex import *
from match_result import MatchResult
from src.primer_data_dto import PrimerDataDTO, get_primer_dto_from_args

# optional: improve offsets to be more accurate
# optional: create options.ini creation + parameterization?

# todo: primers and offsets as files --> check all reads for all pairs with respective offset
# todo: from two lines to biopython parser (since read can be multiline with ~60 char per line)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process sequence alignment parameters.")

    parser.add_argument("--search_area", type=float, default=0.2, help="This value will determine, "
                                                        "how much extra area the smith waterman algorithm will search, "
                                                        "if the other primer has already been found with enough certainty (set by '--sw_cutoff').")
    parser.add_argument("--sw_score_cutoff", type=float, default=0.8, help="Smith-Waterman score cutoff (default: 0.8)")

    parser.add_argument("--primer_information", type=str, default="./data/primer-information.csv",
                        help="CSV list of forward and reverse primer sequence, as well as the expected distance inbetween.")

    parser.add_argument("--input_file_path", type=str, default="./data/DB.COX1.fna", help="Path to input sequence file")
    parser.add_argument("--output_file_path", type=str, default="./data/primer-finder-13-02.csv",
                        help="Path to output results file")

    return parser.parse_args()


def compute_arguments():
    args = parse_arguments()
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

    print(file)
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


def write_output_to_file(output_path, read_metadata, forward_match, backward_match, read):
    with lock, open(output_path, 'a') as out_file:
        out_file.write(read_metadata.replace('|', ';').replace(',', ';').strip() +
                       f"{forward_match.score};{forward_match.read};{forward_match.start_index};" +
                       f"{backward_match.score};{backward_match.read};{backward_match.start_index};{read}\n")


def compute_regex_match(primer, primer_regex, read):
    score = 0
    read_match = ''
    index, end_index = find_exact_match(primer_regex, read)
    if index != -1:
        score = len(primer) * substitution_function('A', 'A')
        read_match = read[index:end_index]
    return MatchResult(score, read_match, index, end_index)


def compute_smith_waterman(primer, read, skip, skip3, substitution):
    score, _, read_match, index = smith_waterman(primer, read, skip, skip3, substitution)
    return MatchResult(score, read_match, index, index + len(read))


### the main processing function

def process_pair(primer_data: PrimerDataDTO, pair):
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
        f_match = compute_smith_waterman(primer_data.f_primer, read[f_search_interval[0]:f_search_interval[1]], -2, -2, substitution_function)
        f_match.start_index += f_search_interval[0]
        f_match.end_index += f_search_interval[0]

        score_threshold = len(primer_data.f_primer) * substitution_function('A', 'A') * primer_data.sw_score_cutoff

        if (b_match.start_index == -1) and (f_match.score > score_threshold):
            b_search_interval = (f_match.end_index + distance - offset, f_match.end_index + distance + len(primer_data.b_primer) + offset)

    if b_match.start_index == -1:
        b_match = compute_smith_waterman(primer_data.b_primer, read[b_search_interval[0]:b_search_interval[1]], -2, -2, substitution_function)
        b_match.start_index += b_search_interval[0]
        b_match.end_index += b_search_interval[0]

    write_output_to_file(primer_data.output_file_path, read_metadata, f_match, b_match, read)


### main script

currentRead = ""

if __name__ == "__main__":
    # todo possibly try chunking
    args = compute_arguments()
    for i, primer_pair in enumerate(args.primer_data):
        primer_data = get_primer_dto_from_args(args, i)
        print(primer_data)
        # todo separate into different files?
        with open(primer_data.output_file_path, 'w') as output_file:
            output_file.write(
                "BOLD ID;Read ID;Country;Phylum;Class;Order;Family;Genus;Species;f_score;f_match;f_index;b_score;b_match;b_index;read\n"
            )

        pairs = read_pairs(primer_data.input_file_path)
        lock = Lock()
        worker = partial(process_pair, primer_data)
        with Pool(initializer=init, initargs=(lock,)) as pool:
            pool.map(worker, pairs)

    print(f"Output has been written to {primer_data.output_file_path}")
