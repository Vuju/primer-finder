import sys
import time
from multiprocessing import Pool, Lock

from smith_waterman import smith_waterman
from primer_finder_regex import *


# Task 3: Implement additional gzip support


### Parameters

proposed_offset = 200
proposed_end_offset = 280

sw_score_cutoff = 0.8

forward_primer = "GGDACWGGWTGAACWGTWTAYCCHCC"
backward_primer = "CCWGTWYTAGCHGGDGCWATYAC"
input_file_path = './data/DB.COX1.fna'
output_file_path = './data/primer-finder-fixed2.csv'


currentRead = ""
forward_primer_regex = regex_builder(forward_primer)
backward_primer_regex = regex_builder(backward_primer)


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



### the main processing function

def process_pair(pair):
    read_metadata, line = pair
    read = line.strip()
    f_search_interval, b_search_interval = (0, len(read)), (0, len(read))
    f_score = 0
    b_score = 0
    f_read = ''
    b_read = ''

    # first check for exact matches
    f_index, f_end_index = find_exact_match(forward_primer_regex, read)
    if f_index != -1:
        f_score = len(forward_primer) * substitution_function('A', 'A')
        f_read = read[f_index:f_end_index]
        b_search_interval = (f_index + proposed_offset, f_index + proposed_end_offset)

    b_index, b_end_index = find_exact_match(backward_primer_regex, read[b_search_interval[0]:b_search_interval[1]])
    if b_index != -1:
        b_index += b_search_interval[0]
        b_end_index += b_search_interval[0]
        b_score = len(backward_primer) * substitution_function('A', 'A')
        b_read = read[b_index:b_end_index]
        f_search_interval = (b_end_index - proposed_end_offset, b_end_index - proposed_offset)

    # for each missing exact match, try smith waterman:
    if f_index == -1:
        f_score, f_primer, f_read, f_index = smith_waterman(forward_primer,
                                                            read[f_search_interval[0]:f_search_interval[1]], -2, -2,
                                                            substitution_function)
        f_index += f_search_interval[0]
        if (b_index == -1) and (f_score > (len(forward_primer) * substitution_function('A', 'A'))):
            b_search_interval = (f_index + len(f_read) + proposed_offset, f_index + len(f_read) + proposed_end_offset)

    if b_index == -1:
        b_score, b_primer, b_read, b_index = smith_waterman(backward_primer,
                                                            read[b_search_interval[0]:b_search_interval[1]], -2, -2,
                                                            substitution_function)
        b_index += b_search_interval[0]

    with lock, open(output_file_path, 'a') as out_file:
        out_file.write(read_metadata.replace('|', ';').replace(',', ';').strip() +
                       f"{f_score};{f_read};{f_index};{b_score};{b_read};{b_index};{read}\n")



### main script

if __name__ == "__main__":
    if len(sys.argv) == 5:
        forward_primer = sys.argv[1]
        backward_primer = sys.argv[2]
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
