import ast
import logging
import pandas as pd
import pyhmmer

from Bio.Seq import Seq

logger = logging.getLogger(__name__)

def list_possible_orf(sequence):
    dna = Seq(sequence)
    orf_list = []
    for frame in range(3):
        framed_seq = dna[frame:]

        # adding trailing N, because Bio is pissed if I don't
        framed_seq = add_trailing_n(framed_seq)

        # Translate using invertebrate mitochondrial code
        protein = framed_seq.translate(table=5)
        if '*' in protein:
            orf_list.append(frame)
    return orf_list

def solve_orfs_for_df(df: pd.DataFrame, threshold = 4, upper_threshold = 50):

    remaining_results = df[~df['possible_orfs'].str.contains(r'\[\]')]
    solved_results = remaining_results[~remaining_results['possible_orfs'].str.contains(',')]
    remaining_results = remaining_results[~(remaining_results['Read ID'].isin(solved_results['Read ID']))]
    solved_results["ORF"] = solved_results.apply(lambda x: x["possible_orfs"][1], axis=1)


    last_solved = pd.DataFrame
    failed = 0
    taxonomic_levels = ['Family', 'Order', 'Class']

    while remaining_results.size > 0:
        current_entry = remaining_results.iloc[0]

        # Try to find enough reference matches at each taxonomic level, from specific to general
        for level in taxonomic_levels:
            level_value = current_entry[level]
            comp_group = solved_results[solved_results[level] == level_value]
            group_size = len(comp_group)

            if group_size >= threshold:
                comp_group = comp_group.sample(min(upper_threshold, group_size))
                related_entries = remaining_results[remaining_results[level] == level_value]
                solved = decide_orfs(comp_group, related_entries)
                solved_results = pd.concat([solved_results, solved], ignore_index=True)
                remaining_results = remaining_results[~(remaining_results[level] == level_value)]
                last_solved = solved
                break
            else:
                logger.log(f"{group_size} entries of {level} {level_value} is too small.")

                # If we've tried all levels and none are big enough
                if level == taxonomic_levels[-1]:
                    remaining_results = remaining_results[~(remaining_results['Family'] == current_entry['Family'])]
                    failed += len(related_entries) if 'related_entries' in locals() else 1


    print(f"A total of {failed} entries were impossible to match.")


## helper functions

def decide_orfs(referenceEntries: pd.DataFrame, questionableEntries: pd.DataFrame):
    alphabet = pyhmmer.easel.Alphabet.amino()

    referenceSequences = []
    referenceEntries.apply(lambda x: referenceSequences.append(build_seq_from_pandas_entry(x)), axis=1)
    questionableSequences = []
    questionableEntries.apply(lambda x: questionableSequences.append(process_ambiguous_orf(x)), axis=1)


    longest_ref_seq = max(len(seq) for seq in referenceSequences)
    longest_que_seq = max(len(seq) for seq in questionableSequences)
    longest_seq_len = max(longest_ref_seq, longest_que_seq)

    referenceSequences = pad_sequences(referenceSequences, minimum=longest_seq_len)
    questionableSequences = pad_sequences(questionableSequences, minimum=longest_seq_len)

    dig_questionable_sequences = [seq.digitize(alphabet=alphabet) for seq in questionableSequences]

    msa = pyhmmer.easel.TextMSA(name=b"myMSA", sequences=referenceSequences)
    background = pyhmmer.plan7.Background(alphabet)
    builder = pyhmmer.plan7.Builder(alphabet)
    digital_msa = msa.digitize(alphabet=alphabet)
    hmm, profile, optimized_profile = builder.build_msa(digital_msa, background)

    pipeline = pyhmmer.plan7.Pipeline(alphabet, background)
    pipeline.bias_filter = False

    sequenceBlock = pyhmmer.easel.DigitalSequenceBlock(alphabet=alphabet, iterable=dig_questionable_sequences)
    hits = pipeline.search_hmm(query=hmm, sequences=sequenceBlock)

    modified_entries = pd.DataFrame(columns=questionableEntries.columns)
    questionableEntries.loc[:, 'ORF'] = ''
    for hit in hits:
        [read_id, correct_orf] = hit.name.decode().split("_")
        questionableEntries.loc[questionableEntries['Read ID'] == read_id, 'ORF'] = correct_orf
        modified_entries = pd.concat([modified_entries, questionableEntries.loc[questionableEntries['Read ID'] == read_id].copy()], ignore_index=True)

        # print(f"changed {read_id} from {prev} to {[correct_orf]}")
    print(f"{len(modified_entries)}({len(hits)} hits) of {len(questionableEntries)} resolved and returned.")
    return modified_entries

def build_seq_from_pandas_entry(entry):
    dna = Seq(entry["read"])
    frame = int(entry["ORF"])
    framed_region = dna[frame:]
    framed_region = add_trailing_n(framed_region)
    protein = framed_region.translate(table=5)
    if len(str(protein)) < 2:
        print(f"Skipping too short sequence for {entry["Read ID"]}: {str(protein)}")
        return None
    text_seq = pyhmmer.easel.TextSequence(name=entry["Read ID"].encode(), sequence=(str(protein)))

    return text_seq

def process_ambiguous_orf(entry):
    possible_orfs = ast.literal_eval(entry["possible_orfs"])
    for possible_orf in possible_orfs:
        dna = Seq(entry["read"])
        dna = add_trailing_n(dna)
        protein = dna.translate(table=5)

        text_seq = pyhmmer.easel.TextSequence(name=(entry["Read ID"].encode() + b"_" + str(possible_orf).encode()),
                                             sequence=str(protein))
        return text_seq

def add_trailing_n(sequence):
    remainder = len(sequence) % 3
    if remainder > 0:
        sequence += ('N' * (3 - remainder))
    return sequence

def pad_sequences(sequences, pad_char='X', minimum=0):
    max_length = max(len(seq) for seq in sequences)
    max_length= max(max_length, minimum)

    padded_sequences = [
        pyhmmer.easel.TextSequence(name=seq.name, sequence=(seq.sequence + pad_char * (max_length - len(seq))))
        for seq in sequences
    ]

    return padded_sequences
