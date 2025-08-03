# Primer Finder

A Python tool for efficient primer sequence identification in DNA reads using both regex pattern matching and Smith-Waterman alignment algorithms, with advanced ORF (Open Reading Frame) detection capabilities.

## Overview

Primer Finder locates forward and reverse primer sequences within DNA reads. It first attempts to find exact matches using regex patterns, then falls back to Smith-Waterman algorithm for more flexible matching when necessary. The tool is optimized with multiprocessing for handling large datasets and includes sophisticated ORF detection using Hidden Markov Models.

## Features

- Fast exact matching using regex patterns
- Sensitive matching with Smith-Waterman alignment algorithm
- Expected distance constraints between primer pairs
- Advanced ORF detection using Hidden Markov Models
- Taxonomic-based ORF resolution
- Configuration-based architecture with YAML support
- Environment variable overrides for configuration
- Multiprocessing support for improved performance
- Comprehensive logging system
- Robust error handling throughout the codebase

## Installation

Clone the repository:

```bash
git clone https://github.com/vuju/primer-finder.git
cd primer-finder
```

Install dependencies:

```bash
pip install -r requirements.txt
```

Alternatively, you can use the environment.yaml file to create a conda environment:

```bash
conda env create -f environment.yml
conda activate primer-finder
```

To use the ORF finder module, you will also need **"muscle"** for Multiple Sequence Alignment (MSA). Please use a version of 5.X with parameter names "-align" and "-output".

## Usage

### Run as script
Run the tool with default configuration:

```bash
python -m primer_finder.cli
```
### Install
Or install it like so:
```bash
pip install -e .

# now you can run it everywhere more easily with
primer-finder
```

You can specify which components to run:

```bash
# Run only primer finding
primer-finder --find-primers
# Or, as an example without installation:
python -m primer_finder.cli --find-primers

# Run only ORF detection
primer-finder --find-orfs

# Run both (default)
primer-finder --find-primers --find-orfs
```

Specify input and output files:

```bash
primer-finder --input path/to/input.fna --output path/to/output.csv    # or
primer-finder --find-primers --input path/to/eyeBOLD.db --table_name specimen
```

## Configuration

Primer Finder uses a configuration system based on YAML files. The default configuration is located at `primer_finder/config/default_config.yaml`. You can create a custom configuration file based on this template.

### Configuration Sections

- **paths**: File paths for input, muscle executable, and log file.
- **database**: Database connection settings including table name, column names, and batch size
- **logging**: Logging settings
- **features**: Feature toggles for primer finder and ORF finder
- **algorithm**: Algorithm parameters including:
  - search_area: Limit for searching second primer when first is found
  - gap_penalty: Penalty for single nucleotide gaps
  - triplet_gap_penalty: Penalty for triplet gaps (preserves reading frame)
  - end_of_read_bonus: Bonus for partial matches at read ends
  - orf_matching thresholds: Controls sequence comparison limits
- **parallelization**: Settings for parallel processing
- **query_parameters**: Primer sequences, cutoff values, expected distances, and taxonomic filters

### Environment Variable Overrides

You can override configuration values using environment variables with the prefix `PRIMER_FINDER_`. For example:

```bash
# Override input file path
export PRIMER_FINDER_PATHS__INPUT_FILE="./data/my_eyeBOLD.db"

# Override number of threads
export PRIMER_FINDER_PARALLELIZATION__NUM_THREADS=4
```

## Input Format

### Database Connection

The tool is designed to work with the eyeBOLD database created by https://github.com/TRojaner2013/eyeBOLD. It requires the specimen table from this database, which contains DNA sequences and associated metadata.

The database connection parameters can be configured in the YAML configuration file:
- table_name: The name of the table containing sequences (default: "specimen")
- id_column_name: The column containing sequence identifiers (default: "specimen_id")
- sequence_column_name: The column containing DNA sequences (default: "nuc_san")
- database_batch_size: Number of records to process in each batch (default: 50000)

## Output

When using eyeBOLD, the database is extended with:
- A "primer_matches" table
  - This contains all the data on forward and reverse primer matches
- A "primer_pairs" table, containing:
  - The inter primer sequence
  - The most likely ORF
  - The amino-acid sequence translated from the most likely ORF

## Architecture

### Primer Finder Module

The primer finder module locates primer sequences within DNA reads:
1. Attempts exact matching using regex patterns
2. Falls back to Smith-Waterman algorithm for approximate matching
3. Uses expected distances between primer pairs to optimize search space
4. Extracts the inter-primer region for ORF analysis

### ORF Finder Module

The ORF finder module identifies and resolves ambiguous open reading frames:
1. Identifies all possible ORFs in the inter-primer region
2. Groups sequences by taxonomic levels (Species, Genus, Family, Order, Class)
3. Builds Hidden Markov Models (HMMs) from solved sequences
4. Queries ambiguous sequences against the HMMs to determine the most likely ORF
5. Uses MUSCLE for multiple sequence alignment

### Connectors

The system uses a connector architecture to handle different input/output sources:
- eyeBOLD_connector: Handles usage of eyeBOLD.
- Other connectors can be implemented by extending the base Connector class.

## Performance Optimization

- Multiprocessing with process pools
- Configurable chunk size for parallel processing
- Search space reduction based on expected primer distances
- Regex-first approach before using more expensive alignment algorithms
- Taxonomic-based grouping for efficient ORF resolution
- Enhanced Smith-Waterman algorithm with triplet gap handling

## License

MIT-Licence
