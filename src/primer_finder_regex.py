import re

def regex_builder(initialString):
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


def find_exact_match(regex, read):
    match = re.search(regex, read)
    if match is not None:
        start_index, end_index = match.span()
        return start_index, end_index
    else:
        return -1, -1


