import re

def regex_builder(initialString):
    regexString = ""
    for char in initialString:
        match char:
            # (see reference table: https://en.wikipedia.org/wiki/Nucleic_acid_notation)
            case 'A':
                regexString += "A"
            case 'C':
                regexString += "C"
            case 'G':
                regexString += "G"
            case 'T':
                regexString += "T"
            case 'U':
                regexString += "T"
            case 'W':
                regexString += "[AT]"
            case 'S':
                regexString += "[CG]"
            case 'M':
                regexString += "[AC]"
            case 'K':
                regexString += "[GT]"
            case 'R':
                regexString += "[AG]"
            case 'Y':
                regexString += "[CT]"
            case 'B':
                regexString += "[CGT]"
            case 'D':
                regexString += "[AGT]"
            case 'H':
                regexString += "[ACT]"
            case 'V':
                regexString += "[ACG]"
            case 'N':
                regexString += "[ACGT]"
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


