"""Command-line interface for primer-finder."""

import argparse
import logging
import sys

from primer_finder.config.config_loader import get_config_loader
from primer_finder.connectors.factory import get_connector

config, is_using_default_config = get_config_loader().get_cli_config()
log_level = config["logging"]["level"]
enable_primer_finder = config["features"]["enable_primer_finder"]
enable_orf_decider = config["features"]["enable_orf_finder"]
log_path = config["paths"]["log_file"]
input_file_path = config["paths"]["input_file"]
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
        '-c', '--config',
        type=str,
        help='Path to a custom configuration file.'
    )
    parser.add_argument(
        "--input", "-i", 
        default=input_file_path,
        help=f"Input file path (default: {input_file_path})"
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


def run_primer_finder(connector, logger):
    """Run the primer finder process."""
    logger.info(f"Starting primer-finding process.")
    primer_finder = PrimerFinder(connector=connector)
    primer_finder.find_all_primers()
    logger.info(f"Primer Finder process has finished.")


def run_orf_finder(connector, logger):
    """Run the ORF finder process."""
    logger.info(f"Starting ORF-matching process.")
    orf_decider = OrfDecider(connector=connector)
    orf_decider.solve_all_orfs()
    logger.info(f"ORF matching completed.")


def main():
    """Main entry point for the CLI."""
    args = parse_args()
    logger = setup_logging(args.log, args.log_level)
    if is_using_default_config:
        logger.info("Using default configuration, as no other configuration was provided.")

    if not args.find_primers and not args.find_orfs:
        args.find_primers = True
        args.find_orfs = True
    if not args.input:
        logger.warning("No input file provided.")
        args.input = input("Please enter the input file path: ").strip()
    connector_args = {
        "db_table_name": args.table_name,
    }
    connector = get_connector(input_file_path=args.input, connector_args=connector_args)

    if args.find_primers:
        run_primer_finder(connector, logger)

    if args.find_orfs:
        run_orf_finder(connector, logger)

    return 0


if __name__ == "__main__":
    sys.exit(main())
