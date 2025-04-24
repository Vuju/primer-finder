from primer_finder.config import ConfigurationError
from primer_finder.matching.connectors.db_connector import DbConnector
from primer_finder.matching.connectors.file_connector import FileConnector


def get_connector(input_file_path: str, connector_args):
    """
    Gets the correct type of connector, depending on the input file type.
    :param input_file_path: Path to the input file.
    :param connector_args: Contains additional information, such as output file path or database table name, depending on the input type.
    :return:
    """
    match input_file_path.split(".")[-1]:
        case "fna", "fasta", "gz":
            if "output_file" in connector_args and connector_args["output_file"]:
                return FileConnector(input_file_path, connector_args["output_file"])
            else:
                return FileConnector(input_file_path, "primer_finder_output.csv")
        case "db":
            return DbConnector(input_file_path, connector_args["db_table_name"])
        case unknown:
            raise ConfigurationError(f"Unknown input type: {unknown}.\n"
                                     f"Currently supported: fna, fasta, gz, db")
