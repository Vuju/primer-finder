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
- *NEW* ORF finder Module: will figure out the best orf by querying possible orfs against a HMM

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
To use the orf-finder module, you will also need **"muscle"** for MSA. Please use a version of 5.X with parameter names "-align" and "-output".

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
| Argument | Default | Description |
|----------|---------|-------------|
| `--primer_finder` | True | Flag as false to disable the primer-searching algorithm. |
| `--orf_matching` | True | Flag as false to disable the orf-decision algorithm. |
| `--search_area` | 0.2 | Determines how much extra area the algorithm will search if the other primer has already been found with enough certainty (set by '--sw_cutoff'). |
| `--sw_score_cutoff` | 0.8 | Smith-Waterman score cutoff for accepting a match. |
| `--primer_information` | ./data/primer-information.csv | CSV list of forward and reverse primer sequences, with expected distances. |
| `--muscle_path` | /mnt/c/Users/Me/bin/muscle | Path to the muscle binary/executable. Runs with version 5.3, using 'muscle_path -align tmp_in.fasta -out tmp_out.fasta'. |
| `--input_file_path` | ./data/DB.COX1.fna | Path to input sequence file. |
| `--output_file_path` | ./data/primer-finder-result.csv | Path to output results file. |
| `--orf_matching_threshold` | 10 | Minimum number of similar sequences required to match an orf. |
| `--orf_matching_upper_threshold` | 50 | Limit of similar sequences used to match an orf. |
| `--protein_translation_table` | 5 | Translation table for Bio.Seq translate(). This is used in orf_finder. |
| `--num_threads` | None | Number of threads to use for processing. |
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
     * Primer finder module
       * Try exact matching via regex
       * If needed, apply Smith-Waterman alignment
       * Record results
     * orf finder module
       * Find all trivial ORFs
       * Define Reference Groups
         * MSA with muscle
         * build a HMM
         * query non-trivial orfs against a suitable HMM

## Performance Optimization

- Multiprocessing with process pools
- Search space reduction based on expected primer distances
- Regex-first approach before using more expensive alignment algorithms

## License

(I didn't think about licensing yet)

