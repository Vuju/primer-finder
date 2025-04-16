import ast
import logging
import subprocess

import numpy as np
import pandas as pd
import pyhmmer

from Bio.Seq import Seq
from pyhmmer.easel import TextSequence
from tqdm import tqdm

from src.config.constants import PROTEIN_TRANSLATION_TABLE, MUSCLE_PATH, E_VALUE, ORF_MATCHING_LOWER_THRESHOLD, \
    ORF_MATCHING_UPPER_THRESHOLD

logger = logging.getLogger(__name__)

def list_possible_orf(sequence, translation_table):
    dna = Seq(sequence)
    orf_list = []
    for frame in range(3):
        framed_sequence = dna[frame:]
        framed_sequence = _trim_to_triplet(framed_sequence)
        protein = framed_sequence.translate(table=translation_table)
        if '*' not in protein:          # The "*" character denotes a stop command in the amino-acid sequence.
            orf_list.append(frame)
    return orf_list


class OrfFinder:

    progress_bar: tqdm
    def __init__(self,
                 translation_table = PROTEIN_TRANSLATION_TABLE,
                 muscle_path = MUSCLE_PATH,
                 e_value_threshold = E_VALUE,
                 lower_reference_threshold: int = ORF_MATCHING_LOWER_THRESHOLD,
                 upper_reference_threshold: int = ORF_MATCHING_UPPER_THRESHOLD,
    ):
        """
        an instance of the OrfFinder class. Parameters should be set on initialization.

        :param translation_table: Translation table for dna --> amino acid for the Bio.Seq.translate() method.
        :param muscle_path: system path to your muscle binary.
        :param e_value_threshold: Sets the e-value threshold for querying an orf-candidate against the HMM.
        :param lower_reference_threshold: Defines the least number of sequences necessary to build an HMM.
        :param upper_reference_threshold: Defines the upper limit for the number of sequences used to build an HMM.
        """
        self.translation_table = translation_table
        self.muscle_path = muscle_path
        self.e_value_threshold = e_value_threshold
        self.lower_reference_threshold = lower_reference_threshold
        self.upper_reference_threshold = upper_reference_threshold

    def solve_orfs_for_df(self, df: pd.DataFrame):
        """
        Takes a dataframe of sequences and attempts to solve all ambiguous orfs.

        :param df: The dataframe of all sequences.

        :return: all sequences with a definite or decided orf.
        """
        unsolved_results = df[~df['possible_orfs'].str.contains(r'\[\]')]
        solved_results = unsolved_results[~unsolved_results['possible_orfs'].str.contains(',')]
        unsolved_results = unsolved_results[~(unsolved_results['Read ID'].isin(solved_results['Read ID']))]
        solved_results["ORF"] = solved_results.apply(lambda x: x["possible_orfs"][1], axis=1)

        not_enough_references = 0
        failed = 0

        self.progress_bar = tqdm(total=len(unsolved_results))
        taxonomic_levels = ['Species', 'Genus', 'Family', 'Order', 'Class']

        while unsolved_results.size > 0:
            current_entry = unsolved_results.iloc[0]

            # Try to match at each taxonomic level, from specific to general
            for level in taxonomic_levels:
                level_value = current_entry[level]
                comparison_group = solved_results[solved_results[level] == level_value]
                group_size = len(comparison_group)

                if group_size >= self.lower_reference_threshold:
                    comparison_group = comparison_group.sample(min(self.upper_reference_threshold, group_size))
                    related_entries = unsolved_results[unsolved_results[level] == level_value]
                    hmm = self._construct_hmm(comparison_group)
                    solved = self._query_sequences_against_hmm(hmm, related_entries)
                    failed += len(related_entries) - len(solved)
                    solved_results = pd.concat([solved_results, solved], ignore_index=True)

                    unsolved_results = unsolved_results[~(unsolved_results[level] == level_value)]
                    break
                else:
                    # If we've tried all levels and none are big enough
                    if level == taxonomic_levels[-1]:
                        number_of_unsolvable_entries = len(solved_results[solved_results['Species'] == current_entry['Species']])
                        logger.warning(f"Removing Species '{current_entry['Species']}' with {number_of_unsolvable_entries} members to continue.")
                        unsolved_results = unsolved_results[~(unsolved_results['Species'] == current_entry['Species'])]
                        not_enough_references += number_of_unsolvable_entries


        self.progress_bar.close()
        logger.info(f"A total of {not_enough_references} entries did not have enough references to match. {failed} were not matched successfully.")
        return solved_results


    ## helper functions

    def _construct_hmm(self, reference_entries: pd.DataFrame) -> pyhmmer.plan7.HMM:
        """
        Takes a list of related sequences, aligns them with "muscle", and builds a HMM.

        :param reference_entries: Set of biologically related, but solved sequences.

        :return: An HMM constructed with an MSA of the input sequences.
        """
        amino_alphabet = pyhmmer.easel.Alphabet.amino()

        reference_sequences = np.zeros(shape=len(reference_entries), dtype=pyhmmer.easel.TextSequence)
        reference_sequences.fill(pyhmmer.easel.TextSequence("".encode(),sequence=""))
        reference_entries = reference_entries.reset_index(drop=True)
        for i, row in reference_entries.iterrows():
            reference_sequences[i] = self._get_amino_text_sequence_of(row)

        # do MSA with muscle for the references
        tmp_in_file = "tmp_in.fasta"
        tmp_out_file = "tmp_out.fasta"
        with open(tmp_in_file, "wb") as f:
            for sequence in reference_sequences:
                # print(seq.sequence)
                sequence.write(f)
        subprocess.run([f"{self.muscle_path}", "-align", tmp_in_file, "-output", tmp_out_file],
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       text=True)

        # build hmm profile from the result
        msa = pyhmmer.easel.MSAFile(tmp_out_file).read()
        msa.name = "tmpMSA".encode()
        background = pyhmmer.plan7.Background(amino_alphabet)
        builder = pyhmmer.plan7.Builder(amino_alphabet)
        digital_msa = msa.digitize(alphabet=amino_alphabet)
        hmm, _, _ = builder.build_msa(digital_msa, background)

        return hmm

    def _query_sequences_against_hmm(self, hmm, ambiguous_entries):
        """
        Takes a set of sequences with multiple possible orfs and queries them against the hmm.
        Returns the (sub)set of decided sequences with a new column "ORF", which contains the most likely orf.

        :param hmm: A Hidden Markov Model constructed with related, but solved sequences.
        :param ambiguous_entries: Set of biologically related sequences with multiple possible orf.

        :return: A set of sequences that have been decided. Decision stored in (new) column "ORF".
        """
        amino_alphabet = pyhmmer.easel.Alphabet.amino()

        ambiguous_sequences = np.zeros(shape=(len(ambiguous_entries), 3), dtype=pyhmmer.easel.TextSequence)
        ambiguous_sequences.fill(pyhmmer.easel.TextSequence("".encode(), sequence=""))
        ambiguous_entries = ambiguous_entries.reset_index(drop=True)
        for i, row in ambiguous_entries.iterrows():
            ambiguous_sequences[i] = self._get_possible_amino_text_sequences_of(row)


        background = pyhmmer.plan7.Background(amino_alphabet)
        pipeline = pyhmmer.plan7.Pipeline(amino_alphabet, background)
        pipeline.bias_filter = False
        pipeline.E = self.e_value_threshold

        ambiguous_entries.loc[:, 'ORF'] = ''
        modified_entries = pd.DataFrame(columns=ambiguous_entries.columns)
        total_hits = 0

        # query all possible orfs against the hmm for each original entry
        for query in ambiguous_sequences:
            digital_ambiguous_sequences = [seq.digitize(alphabet=amino_alphabet) for seq in query]
            sequenceBlock = pyhmmer.easel.DigitalSequenceBlock(alphabet=amino_alphabet,
                                                               iterable=digital_ambiguous_sequences)
            top_hits = pipeline.search_hmm(query=hmm, sequences=sequenceBlock)
            if len(top_hits.reported) > 0:
                total_hits += len(top_hits.reported)
                top_hit = top_hits[0]
                for hit in top_hits:
                    if top_hit is not None and top_hit.name.decode() != "":
                        if hit.evalue < top_hit.evalue:
                            top_hit = hit

                # store best result
                [read_id, correct_orf] = top_hit.name.decode().split("_")
                ambiguous_entries.loc[ambiguous_entries['Read ID'] == read_id, 'ORF'] = correct_orf
                if len(modified_entries[modified_entries["Read ID"] == read_id]) == 0:
                    modified_entries = pd.concat([modified_entries,
                                                  ambiguous_entries.loc[
                                                      ambiguous_entries['Read ID'] == read_id].copy()],
                                                 ignore_index=True)

            if self.progress_bar is not None:
                self.progress_bar.update(1)

        return modified_entries

    def _get_amino_text_sequence_of(self, entry: pd.Series) -> pyhmmer.easel.TextSequence:
        dna = Seq(entry["read"])
        frame = int(entry["ORF"])
        framed_region = dna[frame:]
        framed_region = _trim_to_triplet(framed_region)
        protein = framed_region.translate(table=self.translation_table)
        text_seq = pyhmmer.easel.TextSequence(name=entry["Read ID"].encode(), sequence=(str(protein)))
        return text_seq

    def _get_possible_amino_text_sequences_of(self, entry: pd.Series) -> [TextSequence]:
        possible_orfs = ast.literal_eval(entry["possible_orfs"])
        seqs = np.zeros(shape=3, dtype=pyhmmer.easel.TextSequence)
        seqs.fill(pyhmmer.easel.TextSequence("".encode(), sequence=""))
        for i, possible_orf in enumerate(possible_orfs):
            dna = Seq(entry["read"])
            framed_region = dna[possible_orf:]
            framed_region = _trim_to_triplet(framed_region)
            protein = framed_region.translate(table=self.translation_table)
            text_seq = pyhmmer.easel.TextSequence(name=(entry["Read ID"].encode() + b"_" + str(possible_orf).encode()),
                                                 sequence=str(protein))
            seqs[i] = text_seq
        return seqs


def _trim_to_triplet(sequence):
    remainder = len(sequence) % 3
    if remainder != 0:
        sequence = sequence[:-remainder]
    return sequence