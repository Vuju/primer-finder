import re

def regex_builder(initialString):
    """
    Converts a sequence in nucleic acid notation into a regular expression for DNA, by replacing e.g. "W" with "[AT]".
    :param initialString: A sequence in nucleic acid notation.
    :return:
    """
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
                regexString += "[TU]"
            case 'U':
                regexString += "[TU]"
            case 'W':
                regexString += "[ATU]"
            case 'S':
                regexString += "[CG]"
            case 'M':
                regexString += "[AC]"
            case 'K':
                regexString += "[GTU]"
            case 'R':
                regexString += "[AG]"
            case 'Y':
                regexString += "[CTU]"
            case 'B':
                regexString += "[CGTSKYU]"
            case 'D':
                regexString += "[AGTWKRU]"
            case 'H':
                regexString += "[ACTWMYU]"
            case 'V':
                regexString += "[ACGSMR]"
            case 'N':
                regexString += "."
            case _:
                print(f"Unable to process primer. Didn't find {char}.")
    return regexString


def find_exact_match(regex, read):
    """
    A capsule function for re.search().span().
    Will return -1, -1 instead of None though.
    :param regex: The regex pattern to use.
    :param read: The string to search in.
    :return:
    """
    match = re.search(regex, read)
    if match is not None:
        start_index, end_index = match.span()
        return start_index, end_index
    else:
        return -1, -1


