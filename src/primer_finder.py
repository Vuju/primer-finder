import sys
import time

from multiprocessing import Pool
from smith_waterman import smith_waterman


# Task 3: Implement additional gzip support

forward_primer = "GGDACWGGWTGAACWGTWTAYCCHCC"
backward_primer = "CCWGTWYTAGCHGGDGCWATYAC"
input_file_path = './data/DB.COX1.fna'
output_file_path = './data/primer-finder.csv'

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
            yield (line1, line2)

def process_pair(pair):
    currentRead, line = pair
    read = line.strip()
    f_score, f_primer, f_read, f_index = smith_waterman(forward_primer, read, -2, substitution_function)
    b_score, b_primer, b_read, b_index = smith_waterman(backward_primer, read, -2, substitution_function)
    with open(output_file_path, 'a') as output_file:
        output_file.write(currentRead.replace('|', ';').replace(',', ';').strip() +
                      f"{f_score};{f_primer};{f_index};{b_score};{b_primer};{b_index}\n")

if __name__ == "__main__":
    if len(sys.argv) == 5:
        forward_primer = sys.argv[1]
        backward_primer = sys.argv[2]
        input_file_path = sys.argv[3]
        output_file_path = sys.argv[4]

    with open(output_file_path, 'w') as output_file:
        output_file.write("BOLD ID;Read ID;Country;Phylum;Class;Order;Family;Genus;Species;f_score;f_match;f_index;b_score;b_match;b_index\n")

    pairs = read_pairs(input_file_path)
    with Pool(4) as pool:
        results = pool.map(process_pair, pairs)

    print(f"Output has been written to {output_file_path}")