import os
import sqlite3
from typing import Generator, Tuple, Any

from src.primer_matching.connectors.connector import Connector


class DbConnector(Connector):

    def __init__(self,db_path: str, table_name: str) -> None:
        self.__init_db_connection(db_path)
        self.table_name = table_name

    def get_number_of_sequences(self) -> int:
        pass

    def read_sequences(self) -> Generator[Tuple[str, str], Any, None]:
        pass

    def write_output(self, read_metadata, forward_match, backward_match, inter_primer_sequence, possible_orf):
        pass

    def __init_db_connection(self, db_path: str):
        if not os.path.exists(db_path):
            raise FileNotFoundError(f"File {db_path} does not exist")
        self.db = sqlite3.connect(db_path)
        self.cursor = self.db.cursor()
        if self.cursor.execute(f"SELECT name FROM sqlite_master WHERE name='{self.table_name}'").fetchone() is None:
            raise KeyError(f"No {self.table_name} table in {db_path}")