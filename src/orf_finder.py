import ast
import logging

import numpy as np
import pandas as pd
import pyhmmer

from Bio.Seq import Seq
from tqdm import tqdm

logger = logging.getLogger(__name__)

def list_possible_orf(sequence, translation_table):
    dna = Seq(sequence)
    orf_list = []
    for frame in range(3):
        framed_seq = dna[frame:]

        # trim sequqence, because Bio is pissed if I don't
        framed_seq = trim_to_triplet(framed_seq)

        # Translate using (by default) invertebrate mitochondrial code
        protein = framed_seq.translate(table=translation_table)
        if '*' not in protein:
            orf_list.append(frame)
    return orf_list

def solve_orfs_for_df(
        df: pd.DataFrame,
        translation_table,
        e_value = 1000,
        threshold = 10,
        upper_threshold = 50
):

    remaining_results = df[~df['possible_orfs'].str.contains(r'\[\]')]
    solved_results = remaining_results[~remaining_results['possible_orfs'].str.contains(',')]
    remaining_results = remaining_results[~(remaining_results['Read ID'].isin(solved_results['Read ID']))]
    solved_results["ORF"] = solved_results.apply(lambda x: x["possible_orfs"][1], axis=1)

    failed = 0
    pbar = tqdm(total=len(remaining_results))
    taxonomic_levels = ['Species', 'Genus', 'Family', 'Order', 'Class']

    while remaining_results.size > 0:
        current_entry = remaining_results.iloc[0]

        # Try to match at each taxonomic level, from specific to general
        for level in taxonomic_levels:
            level_value = current_entry[level]
            comp_group = solved_results[solved_results[level] == level_value]
            group_size = len(comp_group)

            if group_size >= threshold:
                comp_group = comp_group.sample(min(upper_threshold, group_size))
                related_entries = remaining_results[remaining_results[level] == level_value]
                solved = decide_orfs_here(comp_group, related_entries,translation_table=translation_table, e_value=e_value, pbar=pbar)
                solved_results = pd.concat([solved_results, solved], ignore_index=True)

                ### the following is a significant speedup over the update(1) updater (30% at 7:20 instead of 9:14),
                ### but of course much less responsive (first update at 7:20 instead of ~2:50)
                # pbar.update(len(related_entries))

                remaining_results = remaining_results[~(remaining_results[level] == level_value)]
                break
            else:
                logger.info(f"{group_size} entries of {level} {level_value} is too small.")

                # If we've tried all levels and none are big enough
                if level == taxonomic_levels[-1]:
                    logger.warning(f"Removing Family '{current_entry['Family']}' with {len(solved_results[solved_results['Family'] == current_entry['Family']])} members to continue.")
                    remaining_results = remaining_results[~(remaining_results['Family'] == current_entry['Family'])]
                    failed += 1 if 'related_entries' not in locals() else len(related_entries)


    pbar.close()
    logger.info(f"A total of {failed} entries were impossible to match.")


## helper functions

def decide_orfs_here(
        referenceEntries: pd.DataFrame,
        questionableEntries: pd.DataFrame,
        translation_table,
        e_value = 1000,
        pbar=None
):
    alphabet = pyhmmer.easel.Alphabet.amino()

    referenceSequences = np.zeros(shape=len(referenceEntries), dtype=pyhmmer.easel.TextSequence)
    referenceSequences.fill(pyhmmer.easel.TextSequence("".encode(),sequence=""))
    referenceEntries = referenceEntries.reset_index(drop=True)
    for i, row in referenceEntries.iterrows():
        referenceSequences[i] = build_seq_from_pandas_entry(row, translation_table=translation_table)

    questionableSequences = np.zeros(shape=(len(questionableEntries), 3), dtype=pyhmmer.easel.TextSequence)
    questionableSequences.fill(pyhmmer.easel.TextSequence("".encode(),sequence=""))
    questionableEntries = questionableEntries.reset_index(drop=True)
    for i, row in questionableEntries.iterrows():
        questionableSequences[i] = process_ambiguous_orf(row, translation_table=translation_table)

    que_lengths = np.vectorize(len)(questionableSequences)
    longest_que_seq = np.max(que_lengths)
    ref_lengths = np.vectorize(len)(referenceSequences)
    longest_ref_seq = np.max(ref_lengths)
    longest_seq_len = max(longest_ref_seq, longest_que_seq)

    referenceSequences = pad_sequences(referenceSequences, minimum=longest_seq_len)
    questionableSequences = pad_sequences_2d(questionableSequences, minimum=longest_seq_len)


    msa = pyhmmer.easel.TextMSA(name=b"myMSA", sequences=referenceSequences.tolist())
    background = pyhmmer.plan7.Background(alphabet)
    builder = pyhmmer.plan7.Builder(alphabet)
    digital_msa = msa.digitize(alphabet=alphabet)
    hmm, profile, optimized_profile = builder.build_msa(digital_msa, background)

    pipeline = pyhmmer.plan7.Pipeline(alphabet, background)
    pipeline.bias_filter = False
    pipeline.E = e_value

    questionableEntries.loc[:, 'ORF'] = ''
    modified_entries = pd.DataFrame(columns=questionableEntries.columns)
    total_hits = 0

    for query in questionableSequences:

        dig_questionable_sequences = [seq.digitize(alphabet=alphabet) for seq in query]
        sequenceBlock = pyhmmer.easel.DigitalSequenceBlock(alphabet=alphabet, iterable=dig_questionable_sequences)
        top_hits = pipeline.search_hmm(query=hmm, sequences=sequenceBlock)
        if len(top_hits.reported) > 0:
            total_hits += len(top_hits.reported)
            top_hit = top_hits[0]

            for hit in top_hits:
                if top_hit is not None and top_hit.name.decode() != "":
                    if hit.evalue < top_hit.evalue:
                        top_hit = hit

            [read_id, correct_orf] = top_hit.name.decode().split("_")
            questionableEntries.loc[questionableEntries['Read ID'] == read_id, 'ORF'] = correct_orf
            if len(modified_entries[modified_entries["Read ID"] == read_id]) == 0:
                modified_entries = pd.concat([modified_entries, questionableEntries.loc[questionableEntries['Read ID'] == read_id].copy()], ignore_index=True)
        ### this is the much slower but more responsive updater (see comment further up)
        if pbar is not None:
            pbar.update(1)
    logger.info(f"modified entries: {len(modified_entries)} (of {total_hits} hits) and {len(questionableEntries)} original entries")
    return modified_entries

def build_seq_from_pandas_entry(entry: pd.Series, translation_table):
    dna = Seq(entry["read"])
    frame = int(entry["ORF"])
    framed_region = dna[frame:]
    framed_region = trim_to_triplet(framed_region)
    protein = framed_region.translate(table=translation_table)

    text_seq = pyhmmer.easel.TextSequence(name=entry["Read ID"].encode(), sequence=(str(protein)))

    return text_seq

def process_ambiguous_orf(entry: pd.Series, translation_table):
    possible_orfs = ast.literal_eval(entry["possible_orfs"])
    seqs = np.zeros(shape=3, dtype=pyhmmer.easel.TextSequence)
    seqs.fill(pyhmmer.easel.TextSequence("".encode(), sequence=""))
    for i, possible_orf in enumerate(possible_orfs):
        dna = Seq(entry["read"])
        framed_region = dna[possible_orf:]
        framed_region = trim_to_triplet(framed_region)
        protein = framed_region.translate(table=translation_table)

        text_seq = pyhmmer.easel.TextSequence(name=(entry["Read ID"].encode() + b"_" + str(possible_orf).encode()),
                                             sequence=str(protein))
        seqs[i] = text_seq
    return seqs

def trim_to_triplet(sequence):
    remainder = len(sequence) % 3
    if remainder != 0:
        sequence = sequence[:-remainder]
    return sequence

def pad_sequences(sequences, minimum=0, pad_char='X'):
    max_length = len(max(sequences, key=len))
    max_length = max(max_length, minimum)

    padded_sequences = np.zeros(shape=len(sequences), dtype=pyhmmer.easel.TextSequence)
    for i, seq in enumerate(sequences):
        padded_sequence = seq.sequence + pad_char * (max_length - len(seq))
        padded_sequences[i] = pyhmmer.easel.TextSequence(name=seq.name, sequence=padded_sequence)

    return np.array(padded_sequences, dtype=pyhmmer.easel.TextSequence)

def pad_sequences_2d(sequences, minimum, pad_char='X'):
    padded_matrix = np.zeros(shape=(len(sequences), len(sequences[0])), dtype=pyhmmer.easel.TextSequence)
    for i in range(len(sequences)):
        padded_matrix[i] = pad_sequences(sequences[i], minimum, pad_char)
    return np.array(padded_matrix)