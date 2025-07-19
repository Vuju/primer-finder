from Bio.Seq import Seq

from primer_finder.orf.decider import _trim_to_triplet


def list_possible_orf(sequence, translation_table):
    """
    A function to determine all reading frames of a sequence, which do not contain a stop codon.
    :param sequence: The sequence to analyze.
    :param translation_table: The translation table for DNA to amino acid.
    :return: returns a List containing 0, 1 and/or 2 as elements, indicating the offset from the start of the sequence as open reading frames.
    """
    dna = Seq(sequence)
    orf_list = []
    for frame in range(3):
        framed_sequence = dna[frame:]
        framed_sequence = _trim_to_triplet(framed_sequence)
        protein = framed_sequence.translate(table=translation_table)
        if '*' not in protein: # The "*" character denotes a stop command in the amino-acid sequence.
            if 'X' not in protein : # And "X" is "unknown"
                orf_list.append(frame)
    return orf_list
