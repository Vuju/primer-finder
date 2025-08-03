import logging
import re

logger = logging.getLogger(__name__)

def regex_builder(initialString: str) -> str:
    """
    Converts a sequence in nucleic acid notation into a regular expression for DNA, by replacing e.g. "W" with "[AT]".
    :param initialString: A sequence in nucleic acid notation.
    :return: A regular expression string or empty string if an error occurs.
    """
    if initialString is None:
        logger.error("Cannot build regex from None input")
        return ""
    
    try:
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
                    regexString += "[WATU]"
                case 'S':
                    regexString += "[SCG]"
                case 'M':
                    regexString += "[MAC]"
                case 'K':
                    regexString += "[KGTU]"
                case 'R':
                    regexString += "[RAG]"
                case 'Y':
                    regexString += "[YCTU]"
                case 'B':
                    regexString += "[BCGTSKYU]"
                case 'D':
                    regexString += "[DAGTWKRU]"
                case 'H':
                    regexString += "[HACTWMYU]"
                case 'V':
                    regexString += "[VACGSMR]"
                case 'N':
                    regexString += "."
                case _:
                    logger.error(f"Unable to process primer. Didn't find {char}.")
                    # Continue building the regex with a wildcard for unknown characters
                    regexString += "."
        return regexString
    except Exception as e:
        logger.error(f"Error building regex from input '{initialString}': {str(e)}")
        return ""


def find_exact_match(regex: str,
                     read: str) -> tuple[int, int]:
    """
    A capsule function for re.search().span().
    Will return -1, -1 instead of None though.
    :param regex: The regex pattern to use.
    :param read: The string to search in.
    :return: A tuple of (start_index, end_index) or (-1, -1) if no match or error.
    """
    if regex is None or read is None:
        logger.error(f"Invalid input for find_exact_match: regex={regex}, read={read}")
        return -1, -1
    
    try:
        match = re.search(regex, read)
        if match is not None:
            start_index, end_index = match.span()
            return start_index, end_index
        else:
            return -1, -1
    except re.error as e:
        logger.error(f"Regular expression error: {str(e)} for pattern '{regex}'")
        return -1, -1
    except Exception as e:
        logger.error(f"Unexpected error in find_exact_match: {str(e)}")
        return -1, -1


