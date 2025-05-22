"""Command-line interface for primer-finder."""

import argparse
import logging
import sys

import pandas as pd

from primer_finder.config.config_loader import get_config_loader
from connectors import get_connector

config = get_config_loader().get_config()
log_level = config["logging"]["level"]
enable_primer_finder = config["features"]["enable_primer_finder"]
enable_orf_decider = config["features"]["enable_orf_finder"]
log_path = config["paths"]["log_file"]
input_file_path = config["paths"]["input_file"]
output_file_path = config["paths"]["output_file"] if "output_file" in config["paths"] else None
db_table_name = config["database"]["table_name"] if "table_name" in config["database"] else None

from primer_finder.matching.primer_finder import PrimerFinder
from primer_finder.orf.decider import OrfDecider


def setup_logging(log_path=log_path, log_level=log_level):
    """Set up the logging configuration."""
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
        default=input_file_path,
        help=f"Input file path (default: {input_file_path})"
    )
    parser.add_argument(
        "--output", "-o", 
        default=output_file_path,
        help=f"Output file path (default: {output_file_path})"
    )
    parser.add_argument(
        "--table_name", "-t",
        default=db_table_name,
        help=f"Name of the db table (default: {db_table_name})"
    )
    parser.add_argument(
        "--find-primers", 
        action="store_true", 
        default=enable_primer_finder,
        help="Run primer finder (default: %(default)s)"
    )
    parser.add_argument(
        "--find-orfs", 
        action="store_true", 
        default=enable_orf_decider,
        help="Run ORF finder (default: %(default)s)"
    )
    parser.add_argument(
        "--log", 
        default=log_path,
        help=f"Log file path (default: {log_path})"
    )
    parser.add_argument(
        "--log-level", 
        default=log_level,
        type=int,
        help=f"Logging level (default: {log_level})"
    )

    return parser.parse_args()


def run_primer_finder(input_file, connector_args, logger):
    """Run the primer finder process."""
    logger.info(f"Starting primer-finding process on {input_file}")

    connector = get_connector(input_file_path=input_file, connector_args=connector_args)
    primer_finder = PrimerFinder(connector=connector)
    primer_finder.find_all_primers()
    logger.info(f"Primer Finder process has finished.")


def run_orf_finder(file_path, logger):
    """Run the ORF finder process."""
    logger.info(f"Starting ORF-matching process")
    orf_finder = OrfDecider()
    all_entries = pd.read_csv(file_path, sep=";")
    solved = orf_finder.solve_orfs_for_df(df=all_entries)
    logger.info("Starting write-back")
    solved.to_csv(file_path)
    logger.info(f"ORF matching output has been written to {file_path}")


def main():
    """Main entry point for the CLI."""
    args = parse_args()
    logger = setup_logging(args.log, args.log_level)

    if not args.find_primers and not args.find_orfs:
        args.find_primers = True
        args.find_orfs = True

    if args.find_primers:
        connector_args = {
            "output_file": args.output,
            "db_table_name": args.table_name,
        }
        run_primer_finder(args.input, connector_args, logger)

    if args.find_orfs:
        run_orf_finder(args.output, logger)

    return 0


if __name__ == "__main__":
    sys.exit(main())
