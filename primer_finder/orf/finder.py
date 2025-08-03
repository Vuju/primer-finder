import logging
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError

from primer_finder.orf.decider import _trim_to_triplet

logger = logging.getLogger(__name__)

def list_possible_orf(sequence, translation_table):
    """
    A function to determine all reading frames of a sequence, which do not contain a stop codon.
    :param sequence: The sequence to analyze.
    :param translation_table: The translation table for DNA to amino acid.
    :return: returns a List containing 0, 1 and/or 2 as elements, indicating the offset from the start of the sequence as open reading frames.
    """
    if sequence is None or len(sequence) == 0:
        logger.warning("Empty or None sequence provided to list_possible_orf")
        return []
    
    try:
        # Validate sequence contains only valid nucleotides
        valid_chars = set("ACGTURYKMSWBDHVN-")
        if not all(c.upper() in valid_chars for c in sequence):
            invalid_chars = [c for c in sequence if c.upper() not in valid_chars]
            logger.warning(f"Sequence contains invalid characters: {invalid_chars[:10]}{'...' if len(invalid_chars) > 10 else ''}")
        
        dna = Seq(sequence)
        orf_list = []
        
        for frame in range(3):
            try:
                framed_sequence = dna[frame:]
                framed_sequence = _trim_to_triplet(framed_sequence)
                
                if len(framed_sequence) == 0:
                    continue
                    
                try:
                    protein = framed_sequence.translate(table=translation_table)
                    if '*' not in protein: # The "*" character denotes a stop command in the amino-acid sequence.
                        if 'X' not in protein: # And "X" is "unknown"
                            orf_list.append(frame)
                except TranslationError as e:
                    logger.error(f"Translation error for frame {frame}: {str(e)}")
                except ValueError as e:
                    logger.error(f"Invalid translation table {translation_table}: {str(e)}")
                    
            except Exception as e:
                logger.error(f"Error processing frame {frame}: {str(e)}")
                
        return orf_list
        
    except Exception as e:
        logger.error(f"Unexpected error in list_possible_orf: {str(e)}")
        return []
