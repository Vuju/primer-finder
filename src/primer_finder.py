import sys
import time

from smith_waterman import smith_waterman


# Task 3: Implement additional gzip support

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


if len(sys.argv) != 5:
    forward_primer = "GGDACWGGWTGAACWGTWTAYCCHCC"
    backward_primer = "CCWGTWYTAGCHGGDGCWATYAC"
    input_file_path = './data/DB.COX1.fna'
    output_file_path = './data/primer-finder.csv'
else:
    forward_primer = sys.argv[1]
    backward_primer = sys.argv[2]
    input_file_path = sys.argv[3]
    output_file_path = sys.argv[4]

currentRead = ""

with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:
    count = 1
    start_time = time.time()
    for line_number, line in enumerate(input_file, start=1):
        if line_number % 2 == 1:
            currentRead = line
        else:
            read = line.strip()
            f_score, f_primer, f_read, f_index = smith_waterman(forward_primer, read, -2, substitution_function)
            b_score, b_primer, b_read, b_index = smith_waterman(backward_primer, read, -2, substitution_function)
            output_file.write(currentRead.replace('|', ';').replace(',', ';').strip() +
                              f"{f_score};{f_primer};{f_index};{b_score};{b_primer};{b_index}\n")
        if line_number >= count * 200:
            count += 1
            print(f"{line_number} at time {time.time()-start_time}.\n")

print(f"Output has been written to {output_file_path}")
