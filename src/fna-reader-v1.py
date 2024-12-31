# imports

input_file_path = './data/DB.COX1.fna'
output_file_path = 'output.txt'

with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:

    # Iterate through each line in the input file
    for line_number, line in enumerate(input_file, start=1):
        # Count the occurrences of 'A' (case-sensitive)
        count_a = line.count('A')
        # Write the count to the output file
        output_file.write(f"Line {line_number}: {count_a}\n")

print(f"Counts of 'A' have been written to {output_file_path}")

