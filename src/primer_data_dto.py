from dataclasses import dataclass
from typing import Any


@dataclass
class PrimerDataDTO:
    f_primer: str
    f_primer_regex: str
    b_primer: str
    b_primer_regex: str

    distance: int

    search_area: float
    sw_score_cutoff: float
    sw_gap: int
    sw_gap3: int

    input_file_path: str
    output_file_path: str

    translation_table: Any

def get_primer_dto_from_args(args, index):
    out_path_split = args.output_file_path.rpartition(".")
    out_path = "".join([out_path_split[0], f"-{index}.", out_path_split[2]])

    return PrimerDataDTO(
        f_primer = args.primer_data[index]["f_primer"],
        f_primer_regex=args.primer_data[index]["f_primer_regex"],
        b_primer=args.primer_data[index]["b_primer"],
        b_primer_regex=args.primer_data[index]["b_primer_regex"],
        distance=int(args.primer_data[index]["distance"]),
        search_area=args.search_area,
        sw_score_cutoff=args.sw_score_cutoff,
        input_file_path=args.input_file_path,
        output_file_path=out_path,
        translation_table=args.protein_translation_table,
        sw_gap=args.gap_penalty,
        sw_gap3=args.gap3_penalty,
    )

