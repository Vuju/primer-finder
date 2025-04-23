"""Command-line interface for primer-finder."""

import argparse
import logging
import sys

import pandas as pd

from primer_finder.config.constants import (
    LOG_PATH,
    LOG_LEVEL,
    PRIMER_FINDER,
    ORF_FINDER,
    OUTPUT_FILE_PATH,
    INPUT_FILE_PATH,
)
from primer_finder.matching.connectors.file_connector import FileConnector
from primer_finder.matching.primer_finder import PrimerFinder
from primer_finder.orf.decider import OrfDecider


def setup_logging(log_path=LOG_PATH, log_level=LOG_LEVEL):
    """Set up logging configuration."""
    logging.basicConfig(filename=log_path, level=log_level)
    logging.getLogger().addHandler(logging.StreamHandler())
    return logging.getLogger(__name__)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Find primers in genetic sequences and analyze open reading frames."
    )
    parser.add_argument(
        "--input", "-i", 
        default=INPUT_FILE_PATH,
        help=f"Input file path (default: {INPUT_FILE_PATH})"
    )
    parser.add_argument(
        "--output", "-o", 
        default=OUTPUT_FILE_PATH,
        help=f"Output file path (default: {OUTPUT_FILE_PATH})"
    )
    parser.add_argument(
        "--find-primers", 
        action="store_true", 
        default=PRIMER_FINDER,
        help="Run primer finder (default: %(default)s)"
    )
    parser.add_argument(
        "--find-orfs", 
        action="store_true", 
        default=ORF_FINDER,
        help="Run ORF finder (default: %(default)s)"
    )
    parser.add_argument(
        "--log", 
        default=LOG_PATH,
        help=f"Log file path (default: {LOG_PATH})"
    )
    parser.add_argument(
        "--log-level", 
        default=LOG_LEVEL,
        type=int,
        help=f"Logging level (default: {LOG_LEVEL})"
    )
    
    return parser.parse_args()


def run_primer_finder(input_file, output_file, logger):
    """Run the primer finder process."""
    logger.info(f"Starting primer-finding process on {input_file}")

    
    connector = FileConnector(input_file=input_file, output_file=output_file)
    primer_finder = PrimerFinder(connector=connector)
    primer_finder.find_all_primers()
    logger.info(f"Primer Finder output has been written to {output_file}")


def run_orf_finder(output_file, logger):
    """Run the ORF finder process."""
    logger.info(f"Starting ORF-matching process")
    orf_finder = OrfDecider()
    all_entries = pd.read_csv(output_file, sep=";")
    solved = orf_finder.solve_orfs_for_df(df=all_entries)
    logger.info("Starting write-back")
    solved.to_csv(output_file)
    logger.info(f"ORF matching output has been written to {output_file}")


def main():
    """Main entry point for the CLI."""
    args = parse_args()
    logger = setup_logging(args.log, args.log_level)
    
    if not args.find_primers and not args.find_orfs:
        logger.error("No actions specified. Use --find-primers or --find-orfs")
        return 1
    
    if args.find_primers:
        run_primer_finder(args.input, args.output, logger)
    
    if args.find_orfs:
        run_orf_finder(args.output, logger)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())