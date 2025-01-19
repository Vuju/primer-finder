import sys
from smith_waterman import smith_waterman

# Task 1: Implement write-back
# Task 2 : Implement proper substitution function

def substitution_function():
    pass


if len(sys.argv) != 5:
    forward_primer = "GGDACWGGWTGAACWGTWTAYCCHCC"
    backward_primer = "CCWGTWYTAGCHGGDGCWATYAC"
    input_file_path = './data/DB.COX1.fna'
    output_file_path = './data/double-primer-matching.txt'
else:
    forward_primer = sys.argv[1]
    backward_primer = sys.argv[2]
    input_file_path = sys.argv[3]
    output_file_path = sys.argv[4]


currentRead = ""

with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:
    for line_number, line in enumerate(input_file, start=1):
        if line_number % 2 == 1:
            currentRead = line
        else:
            f_score, f_primer, f_read, f_index = smith_waterman(forward_primer, line, -2, substitution_function)
            b_score, b_primer, b_read, b_index = smith_waterman(backward_primer, line, -2, substitution_function)

print(f"Output has been written to {output_file_path}")
