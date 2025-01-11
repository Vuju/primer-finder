import string
import sys
import re

def regexBuilder(initialString):
    regexString = ""
    for char in initialString:
        match char:
            # skipping USMKRBVN since they dont appear in my primers (yet) (https://en.wikipedia.org/wiki/Nucleic_acid_notation)
            case 'A':
                regexString += "A"
            case 'T':
                regexString += "T"
            case 'C':
                regexString += "C"
            case 'G':
                regexString += "G"
            case 'W':
                regexString += "[AT]"
            case 'Y':
                regexString += "[CT]"
            case 'D':
                regexString += "[AGT]"
            case 'H':
                regexString += "[ACT]"
            case _:
                print(f"Unable to process primer. Didn't find {char}.")
    return regexString


if len(sys.argv) != 5:
    forward_primer = "GGDACWGGWTGAACWGTWTAYCCHCC"
    backward_primer = "CCWGTWYTAGCHGGDGCWATYAC"
    input_file_path = './data/DB.COX1.fna'
    output_file_path = './data/output.txt'
else:
    forward_primer = sys.argv[1]
    backward_primer = sys.argv[2]
    input_file_path = sys.argv[3]
    output_file_path = sys.argv[4]

forward_primer_regex = regexBuilder(forward_primer)
print(forward_primer_regex)
backward_primer_regex = regexBuilder(backward_primer)
print(backward_primer_regex)
primers = forward_primer_regex + ".*" + backward_primer_regex

currentRead = ""

with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:
    for line_number, line in enumerate(input_file, start=1):
        if line_number % 2 == 1:
            currentRead = line
        else:
            matchAttempt = re.search(primers, line)
            if matchAttempt is not None:
                output_file.write(currentRead.replace('|', ';').replace(',', ';').strip() +
                                  "".join([str(index) + ';' for index in matchAttempt.span()]) +
                                  "\n")


print(f"Output has been written to {output_file_path}")
