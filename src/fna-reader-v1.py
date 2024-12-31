import sys

if len(sys.argv) != 3:
    input_file_path = './data/DB.COX1.fna'
    output_file_path = './data/output.txt'

else:
    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]


with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:

    for line_number, line in enumerate(input_file, start=1):
        count_a = line.count('A')
        output_file.write(f"Line {line_number}: {count_a}\n")

print(f"Counts of 'A' have been written to {output_file_path}")

