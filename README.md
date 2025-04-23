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
- Handles compressed (gzip) and uncompressed input files
- Comprehensive logging system

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

# Run both
primer-finder --find-primers --find-orfs
```

Specify input and output files:

```bash
primer-finder --find-primers --input path/to/input.fna --output path/to/output.csv
```

## Configuration

Primer Finder uses a configuration system based on YAML files. The default configuration is located at `primer_finder/config/default_config.yaml`.

### Configuration Sections

- **paths**: File paths for input, output, primer information, etc.
- **logging**: Logging settings
- **features**: Feature toggles for primer finder and ORF finder
- **algorithm**: Algorithm parameters
- **parallelization**: Settings for parallel processing

### Environment Variable Overrides

You can override configuration values using environment variables with the prefix `PRIMER_FINDER_`. For example:

```bash
# Override input file path
export PRIMER_FINDER_PATHS__INPUT_FILE="./data/my_sequences.fna"

# Override number of threads
export PRIMER_FINDER_PARALLELIZATION__NUM_THREADS=4
```

## Input Format

### Primer Information CSV

The primer information file should be a CSV with the following format:

```
FORWARD_PRIMER_SEQ,REVERSE_PRIMER_SEQ,EXPECTED_DISTANCE
```

Example:
```
GGTCAACAAATCATAAAGATATTGG,TAAACTTCAGGGTGACCAAAAAATCA,650
CCAGAGATTAGAGGGAACTGGATGA,GGGACGGTAAATCATTCAATATTATC,475
```

### Sequence File

The tool supports FASTA format files and can handle both single-line and multi-line sequence entries. Both compressed (.gz) and uncompressed files are supported.

## Output

The output is a CSV file containing:
- Sequence metadata (BOLD ID, Read ID, taxonomic information)
- Forward primer match details (score, matched sequence, position)
- Reverse primer match details
- Inter-primer sequence
- Possible ORFs
- Single, most likely ORF (when ORF finder is used)

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
- FileConnector: Handles file-based input/output
- Other connectors can be implemented by extending the base Connector class

## Performance Optimization

- Multiprocessing with process pools
- Configurable chunk size for parallel processing
- Search space reduction based on expected primer distances
- Regex-first approach before using more expensive alignment algorithms
- Taxonomic-based grouping for efficient ORF resolution

## License

MIT-Licence
