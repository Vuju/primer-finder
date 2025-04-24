import os
import sqlite3
from typing import Generator, Tuple, Any

from primer_finder.config import get_config_loader
from primer_finder.matching.connectors.base import Connector


class DbConnector(Connector):

    def __init__(self,db_path: str, table_name: str = None) -> None:
        config = get_config_loader().get_config()
        self.table_name = table_name or config["database"]["table_name"]

        # todo: parameterize input column names
        self.number_of_sequences = None
        self.__init_db_connection(db_path)

    def get_number_of_sequences(self) -> int:
        if self.number_of_sequences is not None:
            return self.number_of_sequences

        query = f"""
                SELECT COUNT(*) FROM {self.table_name}
                """
        self.cursor.execute(query)
        self.number_of_sequences = self.cursor.fetchone()[0]
        return self.number_of_sequences

    def read_sequences(self, batch_size=1000) -> Generator[Tuple[str, str], Any, None]:
        query = f"""
                SELECT specimenid, nuc_san FROM {self.table_name}
                WHERE forward_match_score IS NULL
                LIMIT {batch_size}
                """

        offset = 0
        while True:
            pagination_query = f"{query} OFFSET {offset}"
            self.cursor.execute(pagination_query)
            rows = self.cursor.fetchall()
            if not rows:
                break
            for specimen_id, sequence in rows:
                yield specimen_id, sequence
            offset += batch_size

    def write_output(self, _, read_metadata, forward_match, backward_match,
                     inter_primer_sequence, possible_orf):
        read_id = read_metadata.split()[0]

        f"{forward_match.score};{forward_match.read};{forward_match.start_index};"
        f"{backward_match.score};{backward_match.read};{backward_match.start_index};"
        f"{inter_primer_sequence};{str(possible_orf)}\n"
        query = f"""
                UPDATE ?
                SET forward_match_score = ?,
                    forward_match_read = ?,
                    forward_match_start_index = ?,
                    backward_match_score = ?,
                    backward_match_read = ?,
                    backward_match_start_index = ?,
                    inter_primer_sequence = ?,
                    possible_orf = ?
                WHERE "specimenid" = ?
                """
        self.cursor.execute(query, (self.table_name,
                                    forward_match.score,
                                    forward_match.read,
                                    forward_match.start_index,
                                    backward_match.score,
                                    backward_match.read,
                                    backward_match.start_index,
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
            "match_score-{sequence-name}": "FLOAT",
            "primer_start_index": "INTEGER",
            "primer_end_index": "INTEGER",
            #"inter_primer_sequence-{sequence-name}-{sequence-name}": "TEXT",
            "most-likely-orf-{sequence-name}-{sequence-name}": "INTEGER"
        }

        # Add missing columns
        for column, data_type in required_columns.items():
            if column not in existing_columns:
                self.cursor.execute(f"ALTER TABLE {self.table_name} ADD COLUMN {column} {data_type}")

        self.db.commit()