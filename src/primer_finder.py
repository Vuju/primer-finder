import sys
from multiprocessing import Pool, Lock

from smith_waterman import smith_waterman
from primer_finder_regex import *
from match_result import MatchResult

# optional: Implement additional gzip support
# optional: improve offsets to be more accurate


### Parameters

proposed_offset = 200
proposed_end_offset = 280

sw_score_cutoff = 0.8

f_primer = "GGDACWGGWTGAACWGTWTAYCCHCC"
b_primer = "CCWGTWYTAGCHGGDGCWATYAC"
input_file_path = './data/DB.COX1.fna'
output_file_path = './data/primer-finder-fixed3.csv'

currentRead = ""
f_primer_regex = regex_builder(f_primer)
b_primer_regex = regex_builder(b_primer)


### Function definitions

def substitution_function(p, r):
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
    with open(file_path, 'r') as file:
        while True:
            line1 = file.readline()
            line2 = file.readline()
            if not line1 or not line2:
                break
            yield line1, line2


def init(l):
    global lock
    lock = l


def write_output_to_file(read_metadata, forward_match, backward_match, read):
    with lock, open(output_file_path, 'a') as out_file:
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

def process_pair(pair):
    read_metadata, line = pair
    read = line.strip()
    f_search_interval, b_search_interval = (0, len(read)), (0, len(read))

    # first check for exact matches
    f_match = compute_regex_match(f_primer, f_primer_regex, read)
    if f_match.start_index != -1:
        b_search_interval = (f_match.start_index + proposed_offset, f_match.start_index + proposed_end_offset)

    b_match = compute_regex_match(b_primer, b_primer_regex, read[b_search_interval[0]:b_search_interval[1]])
    if b_match.start_index != -1:
        b_match.start_index += b_search_interval[0]
        b_match.end_index += b_search_interval[0]
        f_search_interval = (b_match.end_index - proposed_end_offset, b_match.end_index - proposed_offset)

    # for each missing exact match, try smith waterman:
    if f_match.start_index == -1:
        f_match = compute_smith_waterman(f_primer, read[f_search_interval[0]:f_search_interval[1]], -2, -2, substitution_function)
        f_match.start_index += f_search_interval[0]
        f_match.end_index += f_search_interval[0]

        score_threshold = len(f_primer) * substitution_function('A', 'A') * sw_score_cutoff

        if (b_match.start_index == -1) and (f_match.score > score_threshold):
            b_search_interval = (f_match.start_index + proposed_offset, f_match.start_index + proposed_end_offset)

    if b_match.start_index == -1:
        b_match = compute_smith_waterman(b_primer, read[b_search_interval[0]:b_search_interval[1]], -2, -2, substitution_function)
        b_match.start_index += b_search_interval[0]
        b_match.end_index += b_search_interval[0]

    write_output_to_file(read_metadata, f_match, b_match, read)


### main script

if __name__ == "__main__":
    if len(sys.argv) == 5:
        f_primer = sys.argv[1]
        b_primer = sys.argv[2]
        input_file_path = sys.argv[3]
        output_file_path = sys.argv[4]

    with open(output_file_path, 'w') as output_file:
        output_file.write(
            "BOLD ID;Read ID;Country;Phylum;Class;Order;Family;Genus;Species;f_score;f_match;f_index;b_score;b_match;b_index;read\n"
        )
    pairs = read_pairs(input_file_path)
    lock = Lock()
    with Pool(initializer=init, initargs=(lock,)) as pool:
        pool.map(process_pair, pairs)

    print(f"Output has been written to {output_file_path}")
