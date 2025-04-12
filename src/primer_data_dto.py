from dataclasses import dataclass
from typing import Any


@dataclass
class PrimerDataDTO:
    forward_primer: str
    forward_primer_regex: str
    backward_primer: str
    backward_primer_regex: str

    distance: int

    search_area: float
    smith_waterman_score_cutoff: float

    input_file_path: str
    output_file_path: str

    translation_table: Any

def get_primer_dto_from_args(args, index):
    out_path_split = args.output_file_path.rpartition(".")
    out_path = "".join([out_path_split[0], f"-{index}.", out_path_split[2]])

    return PrimerDataDTO(
        forward_primer= args.primer_data[index]["f_primer"],
        forward_primer_regex=args.primer_data[index]["f_primer_regex"],
        backward_primer=args.primer_data[index]["b_primer"],
        backward_primer_regex=args.primer_data[index]["b_primer_regex"],
        distance=int(args.primer_data[index]["distance"]),
        search_area=args.search_area,
        smith_waterman_score_cutoff=args.smith_waterman_score_cutoff,
        input_file_path=args.input_file_path,
        output_file_path=out_path,
        translation_table=args.protein_translation_table,
    )

