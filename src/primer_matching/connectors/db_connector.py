import os
import sqlite3
from typing import Generator, Tuple, Any

from src.primer_matching.connectors.connector import Connector


class DbConnector(Connector):

    def __init__(self,db_path: str, table_name: str) -> None:
        self.table_name = table_name
        self.number_of_sequences = None
        self.__init_db_connection(db_path)

    def get_number_of_sequences(self) -> int:
        if self.number_of_sequences is not None:
            return self.number_of_sequences

        query = f"SELECT COUNT(*) FROM {self.table_name}"
        self.cursor.execute(query)
        self.number_of_sequences = self.cursor.fetchone()[0]
        return self.number_of_sequences

    def read_sequences(self) -> Generator[Tuple[str, str], Any, None]:
        query = f"SELECT specimenid, nuc_san FROM {self.table_name}"
        self.cursor.execute(query)
        for row in self.cursor:
            specimen_id, sequence = row
            yield specimen_id, sequence

    def write_output(self, read_metadata, forward_match, backward_match,
                     inter_primer_sequence, possible_orf):
        read_id = read_metadata.split()[0]

        query = f"""
                UPDATE ?
                SET forward_match = ?,
                    backward_match = ?,
                    inter_primer_sequence = ?,
                    possible_orf = ?
                WHERE "specimenid" = ?
                """
        self.cursor.execute(query, (self.table_name,
                                    forward_match, backward_match,
                                    inter_primer_sequence,
                                    possible_orf, read_id))
        self.db.commit()

    def __init_db_connection(self, db_path: str):
        if not os.path.exists(db_path):
            raise FileNotFoundError(f"File {db_path} does not exist")
        self.db = sqlite3.connect(db_path)
        self.cursor = self.db.cursor()
        if self.cursor.execute("SELECT name FROM sqlite_master WHERE name=?", (self.table_name,)).fetchone() is None:
            raise KeyError(f"No {self.table_name} table in {db_path}")
        self.__ensure_output_columns_exist()

    def __ensure_output_columns_exist(self):
        # Get existing columns
        self.cursor.execute(f"PRAGMA table_info({self.table_name})")
        existing_columns = {row[1] for row in self.cursor.fetchall()}

        # Define columns we need
        required_columns = {
            "forward_match": "TEXT",
            "backward_match": "TEXT",
            "inter_primer_sequence": "TEXT",
            "possible_orf": "TEXT"
        }

        # Add missing columns
        for column, data_type in required_columns.items():
            if column not in existing_columns:
                self.cursor.execute(f"ALTER TABLE {self.table_name} ADD COLUMN {column} {data_type}")

        self.db.commit()