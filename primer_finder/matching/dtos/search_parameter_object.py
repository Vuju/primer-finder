from dataclasses import dataclass
from typing import Any

from primer_finder.matching.regex import regex_builder


@dataclass
class SearchParameterObject:
    """
    A data class to represent a primer pair.
    """
    forward_primer: str
    forward_primer_regex: str
    forward_cutoff: float

    reverse_primer: str
    reverse_primer_regex: str
    reverse_cutoff: float

    distance: int

    protein_translation_table: Any
    taxonomic_filter_rank: str = None
    taxonomic_filter_name: str = None

    def __init__(
            self,
            forward_primer: str,
            reverse_primer: str,
            distance: int,
            forward_cutoff: float,
            reverse_cutoff: float,
            protein_translation_table: Any,
            forward_primer_regex: str = None,
            reverse_primer_regex: str = None,
            taxonomic_filter_rank: str = None,
            taxonomic_filter_name: str = None,
    ):
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer
        self.distance = distance
        self.forward_cutoff = forward_cutoff
        self.reverse_cutoff = reverse_cutoff
        self.protein_translation_table = protein_translation_table

        self.forward_primer_regex = forward_primer_regex or regex_builder(forward_primer)
        self.reverse_primer_regex = reverse_primer_regex or regex_builder(reverse_primer)

        self.taxonomic_filter_rank = taxonomic_filter_rank
        self.taxonomic_filter_name = taxonomic_filter_name


def primer_info_from_config(config_object: dict):
    """
    Parses the primer information from the config object.
    :param config_object: A search_parameters dictionary.
    :return: An object of PrimerDataDTO, wrapping the search parameters.
    """
    return SearchParameterObject(
        forward_primer=config_object["forward_primer"],
        reverse_primer=config_object["reverse_primer"],
        distance=config_object["distance"],
        forward_cutoff=config_object["forward_cutoff"],
        reverse_cutoff=config_object["reverse_cutoff"],
        protein_translation_table=config_object["protein_translation_table"],
        taxonomic_filter_rank=config_object["taxonomic_filter_rank"],
        taxonomic_filter_name=config_object["taxonomic_filter_name"]
    )





