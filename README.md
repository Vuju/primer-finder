# Primer Finder

A Python tool for efficient primer sequence identification in DNA reads using both regex pattern matching and Smith-Waterman alignment algorithms.

## Overview

Primer Finder locates forward and reverse primer sequences within DNA reads. It first attempts to find exact matches using regex patterns, then falls back to Smith-Waterman algorithm for more flexible matching when necessary. The tool is optimized with multiprocessing for handling large datasets.

## Features

- Fast exact matching using regex patterns
- Sensitive matching with Smith-Waterman alignment algorithm
- Expected distance constraints between primer pairs
- Configurable search parameters
- Multiprocessing support for improved performance
- Handles compressed (gzip) and uncompressed input files

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

## Usage

Run the script with default parameters:

```bash
python primer_finder.py
```

Or customize with command-line arguments:

```bash
python primer_finder.py --input_file_path ./data/my_sequences.fna --primer_information ./data/my_primers.csv
```

### Command-line Arguments

| Argument | Default                         | Description                                                                                   |
|----------|---------------------------------|-----------------------------------------------------------------------------------------------|
| `--search_area` | 0.2                             | Determines how much extra area the algorithm will search if the other primer has already been found |
| `--sw_score_cutoff` | 0.8                             | Smith-Waterman score cutoff for accepting a match                                             |
| `--primer_information` | ./data/primer-information.csv   | CSV list of forward and reverse primer sequences, with expected distances                     |
| `--input_file_path` | ./data/DB.COX1.fna              | Path to input sequence file                                                                   |
| `--output_file_path` | ./data/primer-finder-result.csv | Path to output results file                                                                   |

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

The tool supports FASTA format files and can handle both single-line and multi-line sequence entries.

## Output

The output is a CSV file containing:
- Sequence metadata
- Forward primer match details (score, matched sequence, position)
- Reverse primer match details
- Full read sequence

## Processing pipeline:
  1. Parse input arguments
  2. Load primer information
  3. Process each primer pair against all reads
  4. For each read:
     - Try exact matching via regex
     - If needed, apply Smith-Waterman alignment
     - Record results

## Performance Optimization

- Multiprocessing with process pools
- Search space reduction based on expected primer distances
- Regex-first approach before using more expensive alignment algorithms

## License

(I didn't think about licensing yet)

