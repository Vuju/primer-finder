from Bio.Seq import Seq


def find_orf(sequence):
    possibilities = list_possible_orf(sequence)
    if possibilities:
        return possibilities[0]
    return -1

def list_possible_orf(sequence):
    dna = Seq(sequence)
    orf_list = []
    for frame in range(3):
        framed_seq = dna[frame:]
        # Translate using invertebrate mitochondrial code
        protein = framed_seq.translate(table=5)
        if '*' in protein:
            orf_list.append(frame)
    return orf_list
