import sys

if len(sys.argv) != 4:
    forward_primer = "GGDACWGGWTGAACWGTWTAYCCHCC"
    backward_primer = "CCWGTWYTAGCHGGDGCWATYAC"
    primers = [forward_primer, backward_primer]
    input_file_path = './data/DB.COX1.fna'
    output_file_path = './data/output.txt'

else:
    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]


with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:

    for line_number, line in enumerate(input_file, start=1):
        if line_number % 2 == 0:  # Check only Reads, not infos
            if any(substring in line for substring in primers):
                output_file.write(f"Line {line_number}: {line.strip()}\n")
            else:
                output_file.write(f"Line {line_number}: No Primer.\n")
        else:
            output_file.write(line.strip())

print(f"Counts of 'A' have been written to {output_file_path}")

