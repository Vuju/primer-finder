from dataclasses import dataclass
from typing import Any

from src.primer_finder_regex import regex_builder


@dataclass
class PrimerDataDTO:
    """
    A data class to represent a primer pair.
    """
    forward_primer: str
    forward_primer_regex: str
    backward_primer: str
    backward_primer_regex: str

    distance: int

    def __init__(
            self,
            forward_primer: str,
            backward_primer: str,
            distance: int,
            forward_primer_regex: str = None,
            backward_primer_regex: str = None,
    ):
        self.forward_primer = forward_primer
        self.backward_primer = backward_primer
        self.distance = distance

        self.forward_primer_regex = forward_primer_regex or regex_builder(forward_primer)
        self.backward_primer_regex = backward_primer_regex or regex_builder(backward_primer)


def primer_info_from_string(primer_info_string: str):
    """
    A factory function that returns a PrimerDataDTO.
    :param primer_info_string: A comma separated string of a forward primer, backward primer and their most likely distance.
    :return:
    """
    split = primer_info_string.split(",")
    return PrimerDataDTO(
        forward_primer=split[0].strip(),
        backward_primer=split[1].strip(),
        distance=int(split[2].strip()),
    )





