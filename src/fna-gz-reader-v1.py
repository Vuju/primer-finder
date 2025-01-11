import sys
import gzip

if len(sys.argv) != 3:
    input_file_path = './data/DB.COX1.fna.gz'
    output_file_path = './data/output.txt.gz'

else:
    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]


with gzip.open(input_file_path, 'rt') as input_file, gzip.open(output_file_path, 'wb') as output_file:

    for line_number, line in enumerate(input_file, start=1):
        count_a = line.count('A')
        output_file.write(f"Line {line_number}: {count_a}\n".encode('utf-8'))

print(f"Counts of 'A' have been written to {output_file_path}")

