import ast
import logging
import subprocess

import numpy as np
import pandas as pd
import pyhmmer

from Bio.Seq import Seq
from tqdm import tqdm

from primer_finder.connectors.base import Connector
from primer_finder.config.config_loader import get_config_loader

logger = logging.getLogger(__name__)


class OrfDecider:
    """
    Orf Finder is a class to solve ambiguous open reading frames.
    Given a library of sequences, it can compare ambiguous sequences with certain and solved ones to choose the most likely true orf.
    """
    progress_bar: tqdm
    def __init__(self,
                 connector: Connector,
                 chunk_size: int = None,
                 translation_table = None,
                 muscle_path = None,
                 e_value_threshold = None,
                 lower_reference_threshold: int = None,
                 upper_reference_threshold: int = None,
    ):
        """
        an instance of the OrfFinder class. Parameters should be set on initialization.

        :param connector: A connector class necessary to abstract I/O operations.
        :param chunk_size: Number of sequences to process at once.
        :param translation_table: Translation table for dna --> amino acid for the Bio.Seq.translate() method.
        :param muscle_path: system path to your muscle binary.
        :param e_value_threshold: Sets the e-value threshold for querying an orf-candidate against the HMM.
        :param lower_reference_threshold: Defines the least number of sequences necessary to build an HMM.
        :param upper_reference_threshold: Defines the upper limit for the number of sequences used to build an HMM.
        """
        self.connector = connector

        config = get_config_loader().get_config()
        self.chunk_size = chunk_size or config["parallelization"]["chunk_size"]
        self.muscle_path = muscle_path or config["paths"]["muscle"]
        self.translation_table = translation_table or config["algorithm"]["protein_translation_table"]
        self.e_value_threshold = e_value_threshold or config["algorithm"]["e_value"]
        self.lower_reference_threshold = lower_reference_threshold or config["algorithm"]["orf_matching_lower_threshold"]
        self.upper_reference_threshold = upper_reference_threshold or config["algorithm"]["orf_matching_upper_threshold"]

    def solve_all_orfs(self):
        # chunk-wise iterate over all primer pairs and set trivial ones.
        # connector: read chunk
        # for seq in chunk:
        # check if trivial
        # prepare writeback if trivial
        # connector: writeback chunk
        # connector: while next exists: get the next open sequence
        # connector: search for a sufficient comparison group
        # build hmm
        # connector: fetch batch of related sequences
        # query chunk of sequences against hmm
        # connector: writeback results
        """
        Process and solve all ORFs using the provided connector for data operations.
        """
        # Statistics tracking
        total_processed = 0
        trivially_solved = 0
        hmm_solved = 0
        not_enough_references = 0
        failed = 0

        # Process trivial cases first in chunks
        # chunk_generator = self.connector.read_pairs_chunk(self.chunk_size)
        # for chunk in chunk_generator:
        #     solved = self._process_trivial_orfs(chunk)
        #     trivially_solved += len(solved)
        #     total_processed += len(chunk)
        #     if not solved.empty:
        #         # todo possibly collect for fewer writes
        #         self.connector.write_pair_chunk(solved)
        #
        # logger.info(f"Processed {total_processed} sequences. Trivially solved: {trivially_solved}")

        remaining_count = self.connector.get_remaining_unsolved_count()
        self.progress_bar = tqdm(total=remaining_count)

        taxonomic_levels = ['taxon_species', 'taxon_genus', 'taxon_family', 'taxon_order', 'taxon_class']

        # Process non-trivial cases using HMM
        while True:
            current_entry = self.connector.get_next_unsolved_sequence()
            if current_entry is None:
                break

            solved_this_iteration = False

            for level in taxonomic_levels:
                comparison_group, success = self.connector.fetch_sampled_solved_related_sequences(
                    current_entry,
                    level,
                    self.lower_reference_threshold,
                    self.upper_reference_threshold
                )

                if success:
                    hmm = self._construct_hmm(comparison_group)
                    # todo: chunk process this:
                    related_entries = self.connector.fetch_unsolved_related_sequences(current_entry, level)
                    solved = self._query_sequences_against_hmm(hmm, related_entries)

                    hmm_solved += len(solved)
                    failed += len(related_entries) - len(solved)
                    if self.progress_bar:
                        self.progress_bar.update(len(related_entries))

                    if not solved.empty:
                        self.connector.write_pair_chunk(solved)

                    solved_this_iteration = True
                    break

            if not solved_this_iteration:
                entry_id = current_entry['specimen_id']
                unsolved_related_species = self.connector.fetch_unsolved_related_sequences(current_entry, 'taxon_species')
                unsolvable_count = len(unsolved_related_species)

                logger.warning(f"Removing Species '{entry_id}' and {unsolvable_count} members of the same species to continue.")
                not_enough_references += unsolvable_count

                # Update progress bar
                if self.progress_bar:
                    self.progress_bar.update(unsolvable_count)

        if self.progress_bar:
            self.progress_bar.close()

        logger.info(f"HMM solved: {hmm_solved}")
        logger.info(f"A total of {not_enough_references} entries did not have enough references to match.")
        logger.info(f"{failed} entries were not matched successfully.")

    def solve_orfs_for_df(self, df: pd.DataFrame):
        """
        Takes a dataframe of sequences and attempts to solve all ambiguous orfs.

        :param df: The dataframe of all sequences.

        :return: all sequences with a definite or decided orf.
        """
        unsolved_results = df[~df['possible_orfs'].str.contains(r'\[\]')]
        solved_results = unsolved_results[~unsolved_results['possible_orfs'].str.contains(',')]
        unsolved_results = unsolved_results[~(unsolved_results["specimen_id"].isin(solved_results["specimen_id"]))]
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
                ambiguous_entries.loc[ambiguous_entries["specimen_id"] == int(read_id), 'orf_index'] = int(correct_orf)
                if len(modified_entries[modified_entries["specimen_id"] == int(read_id)]) == 0:
                    modified_entries = pd.concat([modified_entries,
                                                  ambiguous_entries.loc[
                                                      ambiguous_entries["specimen_id"] == int(read_id)].copy()],
                                                 ignore_index=True)

            if self.progress_bar is not None:
                self.progress_bar.update(1)

        return modified_entries

    def _get_amino_text_sequence_of(self, entry: pd.Series) -> pyhmmer.easel.TextSequence:
        aa_sequence = entry["orf_aa"]
        text_seq = pyhmmer.easel.TextSequence(name=str(entry["specimen_id"]).encode(), sequence=aa_sequence)
        return text_seq

    def _get_possible_amino_text_sequences_of(self, entry: pd.Series):
        possible_orfs = _decrypt_po(entry["orf_candidates"])
        seqs = np.zeros(shape=3, dtype=pyhmmer.easel.TextSequence)
        seqs.fill(pyhmmer.easel.TextSequence("".encode(), sequence=""))
        for i, possible_orf in enumerate(possible_orfs):
            aa_sequence = self._dna_to_aa(entry["inter_primer_sequence"], possible_orf)
            text_seq = pyhmmer.easel.TextSequence(name=(str(entry["specimen_id"]).encode() + b"_" + str(possible_orf).encode()),
                                                 sequence=str(aa_sequence))
            seqs[i] = text_seq
        return seqs

    def _process_trivial_orfs(self, chunk):
        chunk["orf_index"] = chunk.apply(
            lambda x: x["orf_candidates"][0] if len(_decrypt_po(x["orf_candidates"])) == 1 else x["orf_index"],
            axis=1)
        chunk["orf_aa"] = chunk.apply(
            lambda x: x["orf_aa"] if pd.isna(x["orf_index"]) else
            self._dna_to_aa(x["inter_primer_sequence"], x["orf_index"]),
            axis=1)
        return chunk

    def _dna_to_aa(self, dna_sequence, offset):
        sequence = Seq(dna_sequence)[offset:]
        trimmed = _trim_to_triplet(sequence)
        translated = trimmed.translate(table=self.translation_table)
        return str(translated)

def _trim_to_triplet(sequence):
    remainder = len(sequence) % 3
    if remainder != 0:
        sequence = sequence[:-remainder]
    return sequence

def _decrypt_po(possible_orf: int) -> list[int]:
    return [i for i in range(possible_orf.bit_length()) if possible_orf & (1 << i)]
