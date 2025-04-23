import logging

## File paths
PRIMER_INFORMATION_PATH = "./primer_finder/config/primer_information.csv"
MUSCLE_PATH = "/mnt/c/Users/JulianVu/bin/muscle"
INPUT_FILE_PATH = "./data/DB.COX1.mini.fna"
OUTPUT_FILE_PATH = "./data/primer-finder-result.csv"

# Logger settings
LOG_PATH = "primer-finder_log.log"
LOG_LEVEL = logging.INFO

# Toggle major functionality on or off.
PRIMER_FINDER = False
ORF_FINDER = True

## Define parameters for primer-finder
# Limit the area to search for a second primer, when the first was found.
# Percentage of the distance between the primers.
SEARCH_AREA = 0.2
# Percentage of maximum score. If a match is better,
# The search area for the other primer will be limited.
SMITH_WATERMAN_SCORE_CUTOFF = 0.8
# Values for Smith-Waterman Scoring:
GAP_PENALTY = -2
TRIPLET_GAP_PENALTY = -2
END_OF_READ_BONUS = 1
# The substitution function should compare two characters and return a score.
# Set to "None" to keep the default function.
CUSTOM_SUBSTITUTION_FUNCTION = None


# Limit the number of sequences that a possible orf will be compared to.
ORF_MATCHING_LOWER_THRESHOLD = 10
ORF_MATCHING_UPPER_THRESHOLD = 50

# Protein translation table for "Bio.Seq.translate". Can be an NCBI index or name.
PROTEIN_TRANSLATION_TABLE = 5

## parallelization settings
# Set to "None" to select the maximum number for your cpu.
CUSTOM_NUM_THREADS = None
CHUNKSIZE = 100

E_VALUE = 1000
