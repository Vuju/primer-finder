import logging
import os
import sqlite3
import random
from typing import Generator, Any

import pandas as pd

from primer_finder.config import get_config_loader
from primer_finder.connectors.base import Connector
from primer_finder.matching.dtos.match_result_dto import MatchResultDTO

logger = logging.getLogger(__name__)

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
        self.cutoff = config["algorithm"]["smith_waterman_score_cutoff"]

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
        result = db.execute(query)
        self.number_of_sequences = result.fetchone()[0]
        db.close()
        return self.number_of_sequences

    def read_sequences(self, forward_primer, backward_primer, batch_size=50000) -> Generator[
        tuple[Any, Any, MatchResultDTO, MatchResultDTO], Any, None]:
        query = f"""
                SELECT 
                    input.{self.input_id_column_name} as specimen_id,
                    input.{self.input_sequence_column_name} as sequence,
                    forward_match.primer_start_index as forward_start,
                    forward_match.primer_end_index as forward_end,
                    forward_match.match_score as forward_score,
                    backward_match.primer_start_index as backward_start,
                    backward_match.primer_end_index as backward_end,
                    backward_match.match_score as backward_score
                FROM 
                    {self.input_table_name} as input
                LEFT JOIN 
                    primer_matches as forward_match
                ON 
                    input.{self.input_id_column_name} = forward_match.specimen_id
                    AND forward_match.primer_sequence = ?
                LEFT JOIN 
                    primer_matches as backward_match
                ON 
                    input.{self.input_id_column_name} = backward_match.specimen_id
                    AND backward_match.primer_sequence = ?
            """

        offset = 0
        db = sqlite3.connect(self.db_path)

        while True:
            # Add pagination to the query
            pagination_query = f"{query} LIMIT {batch_size} OFFSET {offset}"
            # Use pandas to read the results into a DataFrame
            df = pd.read_sql_query(
                sql=pagination_query,
                con=db,
                params=[forward_primer, backward_primer],
            )

            # Break if no more results
            if df.empty:
                break

            # Process each row in the DataFrame
            for _, row in df.iterrows():
                specimen_id = row['specimen_id']
                sequence = row['sequence']

                # Create forward primer match DTO
                forward_match = MatchResultDTO(
                    score=row['forward_score'] if pd.notna(row['forward_score']) else -1,
                    read=sequence,
                    start_index=int(row['forward_start']) if pd.notna(row['forward_start']) else -1,
                    end_index=int(row['forward_end']) if pd.notna(row['forward_end']) else -1,
                    primer_sequence=forward_primer
                )

                # Create backward primer match DTO
                backward_match = MatchResultDTO(
                    score=row['backward_score'] if pd.notna(row['backward_score']) else -1,
                    read=sequence,
                    start_index=int(row['backward_start']) if pd.notna(row['backward_start']) else -1,
                    end_index=int(row['backward_end']) if pd.notna(row['backward_end']) else -1,
                    primer_sequence=backward_primer
                )

                yield specimen_id, sequence, forward_match, backward_match

            # Move to the next batch
            offset += batch_size

        # Close the database connection
        db.close()

    def write_output(self, _, information):
        db = sqlite3.connect(self.db_path)
        db.execute('PRAGMA synchronous = NORMAL')
        db.execute('PRAGMA journal_mode=WAL')
        db.execute('PRAGMA journal_size_limit = 0')
        db.execute('PRAGMA cache_size = -2000')
        db.execute('PRAGMA temp_store = MEMORY')
        cursor = db.cursor()
        try:
            # Start a single transaction for all entries
            cursor.execute("BEGIN TRANSACTION")

            # Prepare lists to collect all values for batch insertion
            primer_matches_data = []
            primer_pairs_data = []

            # Collect data for all entries
            for entry in information:
                # Unpack the tuple elements
                specimen_id, forward_match, backward_match, inter_primer_sequence, possible_orf = entry

                # Add forward primer match data
                primer_matches_data.append((
                    specimen_id,
                    forward_match.primer_sequence,
                    forward_match.start_index,
                    forward_match.end_index,
                    forward_match.score
                ))

                # Add backward primer match data
                primer_matches_data.append((
                    specimen_id,
                    backward_match.primer_sequence,
                    backward_match.start_index,
                    backward_match.end_index,
                    backward_match.score
                ))

            # Execute batch insert for primer matches
            cursor.executemany("""
                INSERT INTO primer_matches
                (specimen_id, primer_sequence, primer_start_index, primer_end_index, match_score)
                VALUES (?, ?, ?, ?, ?)
                    ON CONFLICT(specimen_id, primer_sequence) DO UPDATE 
                    SET match_score = excluded.match_score,
                        primer_start_index = excluded.primer_start_index,
                        primer_end_index = excluded.primer_end_index
            """, primer_matches_data)

            # We need to retrieve the IDs for each pair to link in primer_pairs table
            # This requires us to query back the inserted records
            forward_ids = {}
            backward_ids = {}

            for i, entry in enumerate(information):
                specimen_id, forward_match, backward_match, inter_primer_sequence, possible_orf = entry

                # Get IDs for forward primers
                cursor.execute("""
                               SELECT rowid
                               FROM primer_matches
                               WHERE specimen_id = ?
                                 AND primer_sequence = ?
                               """, (specimen_id, forward_match.primer_sequence))
                forward_match_id = cursor.fetchone()[0]

                # Get IDs for backward primers
                cursor.execute("""
                               SELECT rowid
                               FROM primer_matches
                               WHERE specimen_id = ?
                                 AND primer_sequence = ?
                               """, (specimen_id, backward_match.primer_sequence))
                backward_match_id = cursor.fetchone()[0]

                # Add primer pair data
                primer_pairs_data.append((
                    forward_match_id,
                    backward_match_id,
                    specimen_id,
                    inter_primer_sequence,
                    _encrypt_po(possible_orf),
                    self._set_flags(forward_match, backward_match)
                ))

            # Execute batch insert for primer pairs
            cursor.executemany("""
                INSERT OR REPLACE INTO primer_pairs
                (forward_match_id, reverse_match_id, specimen_id, 
                inter_primer_sequence, orf_candidates, matching_flags)
                VALUES (?, ?, ?, ?, ?, ?)
            """, primer_pairs_data)

            db.commit()
            #wal_size = os.path.getsize("/mnt/z/Uni/Master Thesis/eyeBOLD/eyeBOLD_mini.db-wal")
            #logger.info(f"WAL size: {wal_size / (1024*1024):.2f} MB")

        except Exception as e:
            # Rollback in case of error
            logger.error(f"Failed to write primer pair data: {str(e)}")
            cursor.execute("ROLLBACK")
            raise Exception(f"Failed to write primer pair data: {str(e)}")
        finally:
            cursor.close()
            db.close()

    def _set_flags(self, forward_match, backward_match) -> int:
        f_cutoff = self.cutoff * 2 * len(forward_match.primer_sequence)
        b_cutoff = self.cutoff * 2 * len(backward_match.primer_sequence)
        if forward_match.score < f_cutoff:
            if backward_match.score < b_cutoff:
                return -3
            return -1
        else:
            if backward_match.score < b_cutoff:
                return -2
        return 0

    def __init_db_connection(self):
        if not os.path.exists(self.db_path):
            raise FileNotFoundError(f"File {self.db_path} does not exist")
        db = sqlite3.connect(self.db_path)
        if db.execute("SELECT name FROM sqlite_master WHERE name=?", (self.input_table_name,)).fetchone() is None:
            raise KeyError(f"No {self.input_table_name} table in {self.db_path}")
        db.close()

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
                    FOREIGN KEY (specimen_id) REFERENCES {self.input_table_name}({self.input_id_column_name}),
                    UNIQUE (specimen_id, primer_sequence)
                )
            """)
            db.commit()
            db.close()
            return

        # If the table exists, check for missing columns
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
        db.close()

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
            cursor.execute(f"""
                CREATE TABLE primer_pairs (
                    forward_match_id      INTEGER NOT NULL,
                    reverse_match_id      INTEGER NOT NULL,
                    specimen_id           INTEGER NOT NULL,
                    inter_primer_sequence TEXT,
                    orf_candidates        INTEGER,
                    orf_index             INTEGER,
                    orf_aa                TEXT,
                    matching_flags        INTEGER,
                    PRIMARY KEY (forward_match_id, reverse_match_id),
                    FOREIGN KEY (forward_match_id) REFERENCES primer_matches(match_id),
                    FOREIGN KEY (reverse_match_id) REFERENCES primer_matches(match_id),
                    FOREIGN KEY (specimen_id) REFERENCES {self.input_table_name}({self.input_id_column_name}),
                    CHECK (forward_match_id != reverse_match_id)
                )
            """)
            db.commit()
            db.close()
            return

        # If the table exists, check for missing columns
        cursor.execute("PRAGMA table_info(primer_pairs)")
        existing_columns = {row[1] for row in cursor.fetchall()}

        required_columns = {
            "forward_match_id": "INTEGER NOT NULL",
            "reverse_match_id": "INTEGER NOT NULL",
            "specimen_id": "INTEGER NOT NULL",
            "inter_primer_sequence": "TEXT",
            "orf_candidates": "INTEGER",
            "orf_index": "INTEGER",
            "orf_aa": "TEXT",
            "matching_flags": "INTEGER"
        }
        for column, data_type in required_columns.items():
            if column not in existing_columns:
                cursor.execute(f"ALTER TABLE primer_pairs ADD COLUMN {column} {data_type}")
        cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_orf_index ON primer_pairs(orf_index)")
        db.commit()
        db.close()

    # ------------------------  Connector for ORF Matching ------------------------

    def read_pairs_chunk(self, chunk_size, batch_size = 50000):
        query = f"""
                SELECT forward_match_id, reverse_match_id, specimen_id, inter_primer_sequence, 
                orf_candidates, orf_index, orf_aa, matching_flags
                FROM primer_taxonomic_groups                
                """
        offset = 0
        db = sqlite3.connect(self.db_path)
        while True:
            pagination_query = f"{query} LIMIT {batch_size} OFFSET {offset}"
            offset += batch_size
            df = pd.read_sql_query(
                sql=pagination_query,
                con=db,
            )
            if df.empty:
                break
            yield df
        db.close()

    def write_pair_chunk(self, solved, tmp = True):
        db = sqlite3.connect(self.db_path)
        db.execute("PRAGMA journal_mode = WAL")
        db.execute('PRAGMA temp_store = MEMORY')
        try:
            # Use INSERT OR REPLACE to handle conflicts
            insert_query = f'''
            UPDATE {"primer_taxonomic_groups" if tmp else "primer_pairs"}
            SET orf_index = ?, 
                orf_aa = ?
            WHERE forward_match_id = ? 
              AND reverse_match_id = ?
            '''

            # Convert DataFrame to list of tuples
            update_data = solved[['orf_index', 'orf_aa', 'forward_match_id', 'reverse_match_id']]
            data_tuples = [tuple(row) for row in update_data.itertuples(index=False, name=None)]

            db.execute("BEGIN TRANSACTION")
            db.executemany(insert_query, data_tuples)
            db.commit()
        except Exception as e:
            logger.error(f"Failed to write primer pair data: {str(e)}")
            db.execute("ROLLBACK")
            raise Exception(f"Failed to write primer pair data: {str(e)}")
        finally:
            db.close()

    def get_remaining_unsolved_count(self):
        db = sqlite3.connect(self.db_path)
        result = db.execute(f"""
                SELECT COUNT(*) 
                FROM primer_taxonomic_groups
                WHERE orf_index IS NULL
                """)
        return result.fetchone()[0]

    def get_next_unsolved_sequence(self):
        db = sqlite3.connect(self.db_path)
        try:
            query = """
                    SELECT *
                    FROM primer_taxonomic_groups
                    WHERE orf_index IS NULL
                    LIMIT 1
                    """
            df = pd.read_sql_query(query, db)
            db.close()
            return df
        except Exception as e:
            logger.error(f"Error fetching data: {e}")
            db.close()
            return pd.DataFrame()  # Return empty DataFrame on error

    def fetch_sampled_solved_related_sequences(self, current_entry, level, lower_reference_threshold, upper_reference_threshold, random_seed: int = None):
        if random_seed is not None:
            random.seed(random_seed)

        matching_entries, found_sequences = self._fetch_any_related_sequences(current_entry, level, True)
        if found_sequences:
            filtered_entries  = matching_entries[(matching_entries["orf_index"] >= 0) & (matching_entries["matching_flags"] == 0)]
            if len(filtered_entries) > lower_reference_threshold:
                # Randomly sample up to max_entries from this group
                sample_size = min(upper_reference_threshold, len(filtered_entries))
                selected_entries = filtered_entries.sample(n=sample_size, axis=0, random_state=random_seed)
                return selected_entries, True
        return None, False


    def fetch_unsolved_related_sequences(self, current_entry, level):
        return self._fetch_any_related_sequences(current_entry, level, False)

    def init_temp_pairs_table(self, forward_primer_seq, reverse_primer_seq):
        conn = sqlite3.connect(self.db_path)
        # self.remove_temp_table()
        creation_query = f"""
            CREATE TABLE IF NOT EXISTS primer_taxonomic_groups AS
            SELECT pp.*, s.taxon_genus, s.taxon_species, s.taxon_family, s.taxon_order, s.taxon_class,
                   fm.primer_sequence as forward_seq, rm.primer_sequence as reverse_seq
            FROM primer_pairs pp
            JOIN primer_matches fm ON pp.forward_match_id = fm.match_id
            JOIN primer_matches rm ON pp.reverse_match_id = rm.match_id
            JOIN specimen s ON pp.specimen_id = s.specimenid
            WHERE fm.primer_sequence = ?
            AND rm.primer_sequence = ?;
            """
        logger.info(f"Creating temp pairs table for {forward_primer_seq} and {reverse_primer_seq}. This may take a while.")
        conn.execute(creation_query,(forward_primer_seq, reverse_primer_seq))
        logger.info("Finished creating temp pairs table. Creating indexes.")

        logger.info("Setting up indexes: Orf Index (1/7)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_sequence ON primer_taxonomic_groups(orf_index)")
        logger.info("Setting up indexes: Species (2/7)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_species ON primer_taxonomic_groups(taxon_species)")
        logger.info("Setting up indexes: Genus (3/7)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_genus ON primer_taxonomic_groups(taxon_genus)")
        logger.info("Setting up indexes: Family (4/7)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_family ON primer_taxonomic_groups(taxon_family)")
        logger.info("Setting up indexes: Order (5/7)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_order ON primer_taxonomic_groups(taxon_order)")
        logger.info("Setting up indexes: Class (6/7)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_class ON primer_taxonomic_groups(taxon_class)")
        logger.info("Setting up indexes: Orf Index (7/7)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_match_ids ON primer_taxonomic_groups(forward_match_id, reverse_match_id)")
        logger.info("Finished setting up indexes.")
        conn.close()

    def remove_temp_table(self):
        conn = sqlite3.connect(self.db_path)
        conn.execute("DROP TABLE IF EXISTS primer_taxonomic_groups")
        conn.close()

    def _fetch_any_related_sequences(self, current_entry, level, solved: bool):
        conn = sqlite3.connect(self.db_path)
        try:
            related_query = f"""
                SELECT forward_match_id, reverse_match_id, specimen_id,
                       inter_primer_sequence, orf_candidates, orf_index, 
                       orf_aa, matching_flags
                FROM primer_taxonomic_groups
                WHERE {level} = ?
                AND orf_index {"IS NOT NULL" if solved else "IS NULL"}
                """
            matching_entries = pd.read_sql_query(
                sql=related_query,
                con=conn,
                params=[current_entry[level][0]],
            )
            #logger.info(f"Found {len(matching_entries)} related and {"solved" if solved else "unsolved"} entries in {current_entry[level][0]}")
            return matching_entries, True

        except Exception as e:
            logger.error(f"Error while finding related group for: {e}")
            return None, False
        finally:
            conn.close()

    def primer_pairs_writeback(self):
        conn = sqlite3.connect(self.db_path)
        writeback_query = """
        UPDATE primer_pairs
        SET orf_index = ptg.orf_index,
            orf_aa = ptg.orf_aa
        FROM primer_taxonomic_groups ptg
        WHERE primer_pairs.forward_match_id = ptg.forward_match_id
          AND primer_pairs.reverse_match_id = ptg.reverse_match_id
        """
        conn.execute(writeback_query)
        conn.commit()
        conn.close()
