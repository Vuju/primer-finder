import os
import sqlite3
from typing import Generator, Tuple, Any

from primer_finder.config import get_config_loader
from primer_finder.matching.connectors.base import Connector


def _encrypt_po(possible_orf: list[int]) -> int:
    return sum(1 << orf for orf in possible_orf)

def _decrypt_po(possible_orf: int) -> list[int]:
    return [i for i in range(possible_orf.bit_length()) if possible_orf & (1 << i)]


class DbConnector(Connector):

    def __init__(self,db_path: str, table_name: str = None) -> None:
        config = get_config_loader().get_config()
        self.input_table_name = table_name or config["database"]["input_table_name"]
        self.input_id_column_name = config["database"]["id_column_name"]
        self.input_sequence_column_name = config["database"]["sequence_column_name"]

        self.number_of_sequences = None
        self.db_path = db_path
        self.__init_db_connection()
        self.__ensure_input_table_exists()
        self.__ensure_primer_matches_table_exists()
        self.__ensure_primer_pairs_table_exists()

    def get_number_of_sequences(self) -> int:
        if self.number_of_sequences is not None:
            return self.number_of_sequences

        query = f"""
                SELECT COUNT(*) FROM {self.input_table_name}
                """
        db = sqlite3.connect(self.db_path)
        cursor = db.cursor()
        cursor.execute(query)
        self.number_of_sequences = cursor.fetchone()[0]
        return self.number_of_sequences

    def read_sequences(self, batch_size=1000) -> Generator[Tuple[str, str], Any, None]:
        query = f"""
                SELECT {self.input_id_column_name}, {self.input_sequence_column_name} FROM {self.input_table_name}
                LIMIT {batch_size}
                """

        offset = 0
        db = sqlite3.connect(self.db_path)
        cursor = db.cursor()
        while True:
            pagination_query = f"{query} OFFSET {offset}"
            cursor.execute(pagination_query)
            rows = cursor.fetchall()
            if not rows:
                break
            for specimen_id, sequence in rows:
                yield specimen_id, sequence
            offset += batch_size

    def write_output(self, _, specimen_id, forward_match, backward_match, __, possible_orf):
        db = sqlite3.connect(self.db_path)
        cursor = db.cursor()
        try:
            cursor.execute("BEGIN TRANSACTION")

            cursor.execute("""
                                INSERT INTO primer_matches
                                (specimen_id, primer_sequence, primer_start_index, primer_end_index, match_score)
                                VALUES (?, ?, ?, ?, ?)
                                """, (
                                    specimen_id,
                                    forward_match.primer_sequence,
                                    forward_match.start_index,
                                    forward_match.end_index,
                                    forward_match.score
                                ))
            forward_match_id = cursor.lastrowid
            cursor.execute("""
                                INSERT INTO primer_matches
                                (specimen_id, primer_sequence, primer_start_index, primer_end_index, match_score)
                                VALUES (?, ?, ?, ?, ?)
                                """, (
                                    specimen_id,
                                    backward_match.primer_sequence,
                                    backward_match.start_index,
                                    backward_match.end_index,
                                    backward_match.score
                                ))
            backward_match_id = cursor.lastrowid

            cursor.execute("""
                                INSERT INTO primer_pairs
                                    (forward_match_id, reverse_match_id, orf_candidates)
                                VALUES (?, ?, ?)
                                """, (
                                    forward_match_id,
                                    backward_match_id,
                                    _encrypt_po(possible_orf)
                                ))

            cursor.execute("COMMIT")

        except Exception as e:
            # Rollback in case of error
            cursor.execute("ROLLBACK")
            raise Exception(f"Failed to write primer pair data: {str(e)}")

    def __init_db_connection(self):
        if not os.path.exists(self.db_path):
            raise FileNotFoundError(f"File {self.db_path} does not exist")
        db = sqlite3.connect(self.db_path)
        cursor = db.cursor()
        if cursor.execute("SELECT name FROM sqlite_master WHERE name=?", (self.input_table_name,)).fetchone() is None:
            raise KeyError(f"No {self.input_table_name} table in {self.db_path}")

    def __ensure_input_table_exists(self):
        # todo probably check whether columns exist in input specimen table:
        pass
    def __ensure_primer_matches_table_exists(self):
        db = sqlite3.connect(self.db_path)
        cursor = db.cursor()
        cursor.execute(f"SELECT name FROM sqlite_master WHERE type='table' AND name='primer_matches'")
        table_exists = cursor.fetchone() is not None

        # Create the table if it doesn't exist
        if not table_exists:
            cursor.execute(f"""
                CREATE TABLE primer_matches (
                    match_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    specimen_id INTEGER NOT NULL,
                    primer_sequence TEXT NOT NULL,
                    primer_start_index INTEGER,
                    primer_end_index INTEGER,
                    match_score FLOAT,
                    FOREIGN KEY (specimen_id) REFERENCES {self.input_table_name}({self.input_id_column_name})
                )
            """)
            db.commit()
            return

        # If table exists, check for missing columns
        cursor.execute(f"PRAGMA table_info(primer_matches)")
        existing_columns = {row[1] for row in cursor.fetchall()}

        required_columns = {
            "match_id": "INTEGER PRIMARY KEY AUTOINCREMENT",
            "specimen_id": "INTEGER NOT NULL",
            "primer_sequence": "TEXT NOT NULL",
            "primer_start_index": "INTEGER",
            "primer_end_index": "INTEGER",
            "match_score": "FLOAT"
        }

        for column, data_type in required_columns.items():
            if column not in existing_columns:
                # todo extend for a check for missing primary and foreign key, as those cant be added later
                cursor.execute(f"ALTER TABLE primer_matches ADD COLUMN {column} {data_type}")


        db.commit()

    def __ensure_primer_pairs_table_exists(self):
        """
        Ensures that the primer_pairs table exists in the database.
        This table stores pairs of primer matches that form potential ORFs.
        """
        # Check if the table exists
        db = sqlite3.connect(self.db_path)
        cursor = db.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='primer_pairs'")
        table_exists = cursor.fetchone() is not None

        # Create the table if it doesn't exist
        if not table_exists:
            cursor.execute("""
                CREATE TABLE primer_pairs (
                    forward_match_id INTEGER NOT NULL,
                    reverse_match_id INTEGER NOT NULL,
                    orf_candidates   INTEGER,
                    orf              INTEGER,
                    PRIMARY KEY (forward_match_id, reverse_match_id),
                    FOREIGN KEY (forward_match_id) REFERENCES primer_matches(match_id),
                    FOREIGN KEY (reverse_match_id) REFERENCES primer_matches(match_id),
                    CHECK (forward_match_id != reverse_match_id)
                )
            """)
            db.commit()
            return

        # If table exists, check for missing columns
        cursor.execute("PRAGMA table_info(primer_pairs)")
        existing_columns = {row[1] for row in cursor.fetchall()}

        required_columns = {
            "forward_match_id": "INTEGER NOT NULL",
            "reverse_match_id": "INTEGER NOT NULL",
            "orf_candidates": "INTEGER",
            "orf": "INTEGER"
        }
        for column, data_type in required_columns.items():
            if column not in existing_columns:
                cursor.execute(f"ALTER TABLE primer_pairs ADD COLUMN {column} {data_type}")

        db.commit()