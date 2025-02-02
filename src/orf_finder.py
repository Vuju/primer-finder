from Bio.Seq import Seq

def find_orf(sequence):
    dna = Seq(sequence)
    # Translate using mitochondrial invertebrate code (table 5)
    for frame in range(3):
        framed_seq = dna[frame:]
        # Translate using invertebrate mitochondrial code
        protein = framed_seq.translate(table=5)
        if '*' in protein:
            return frame
    return -1
