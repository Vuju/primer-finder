from primer_finder.config import ConfigurationError
from primer_finder.connectors.eyeBOLD_connector import EyeBOLDConnector


def get_connector(input_file_path: str, connector_args):
    """
    Gets the correct type of connector, depending on the input file type.
    :param input_file_path: Path to the input file.
    :param connector_args: Contains additional information, such as output file path or database table name, depending on the input type.
    :return:
    """
    match input_file_path.split(".")[-1]:
        case "db":
            return EyeBOLDConnector(input_file_path, connector_args["db_table_name"])
        case _:
            raise ConfigurationError(f"Unknown input type: {input_file_path.split(".")[-1]}.\n"
                                     f"Currently supported extension: .db (eyeBOLD: sqlite3)")
