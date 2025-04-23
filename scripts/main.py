import logging

import pandas as pd

from primer_finder.config.constants import LOG_PATH, LOG_LEVEL, PRIMER_FINDER, ORF_FINDER, OUTPUT_FILE_PATH, INPUT_FILE_PATH
from primer_finder.matching.connectors.file_connector import FileConnector
from primer_finder.matching.primer_finder import PrimerFinder
from primer_finder.orf.decider import OrfDecider

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    logging.basicConfig(filename=LOG_PATH, level=LOG_LEVEL)
    logging.getLogger().addHandler(logging.StreamHandler())

    # todo: extract read and write functionality from primer finder

    if PRIMER_FINDER:
        # def connector dynamically
        # def lock
        connector = FileConnector(lock=lock,input_file=INPUT_FILE_PATH,output_file=OUTPUT_FILE_PATH)
        primer_finder = PrimerFinder(connector=connector)
        primer_finder.find_all_primers()
        logger.info(f"Primer Finder output has been written to {OUTPUT_FILE_PATH}")

    # todo make orf finder able to decide for multiple pairs.
        # before, this was handled by having different files. Will be more strange for db.
    if ORF_FINDER:
        logger.info(f"Starting orf-matching process.")
        orf_finder = OrfDecider()
        all_entries = pd.read_csv(OUTPUT_FILE_PATH, sep=";")
        solved = orf_finder.solve_orfs_for_df(df=all_entries)
        logger.info("Starting write-back.")
        solved.to_csv(OUTPUT_FILE_PATH)

        logger.info(f"Orf matching output has been written to {OUTPUT_FILE_PATH}")