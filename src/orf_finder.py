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

        # adding trailing N, because Bio is pissed if I don't
        remainder = len(framed_seq) % 3
        if remainder > 0:
            for i in range(3 - remainder):
                framed_seq += "N"

        # Translate using invertebrate mitochondrial code
        protein = framed_seq.translate(table=5)
        if '*' in protein:
            orf_list.append(frame)
    return orf_list
