import gzip
from multiprocessing import Lock
from typing import TextIO

from primer_finder.matching.connectors.base import Connector


class FileConnector(Connector):

    def __init__(self, input_file, output_file,):
        self.input_file = input_file
        self.output_file = output_file
        self.length = -1

        with open(self.output_file, 'w') as output_file:
            output_file.write(
                "BOLD ID;Read ID;Country;Phylum;Class;Order;Family;Genus;Species;f_score;f_match;f_index;b_score;b_match;b_index;read;possible_orfs\n"
            )

    def get_number_of_sequences(self):
        if self.length != -1:
            return self.length

        count = 0
        file: TextIO
        if self.input_file.endswith('.gz'):
            file = gzip.open(self.input_file, 'rt')
        else:
            file = open(self.input_file, 'r', encoding="UTF-8")
        for line in file.readlines():
            if line.startswith('>'):
                count += 1
        file.close()
        self.length = count
        return self.length

    def read_sequences(self):
        file: TextIO
        if self.input_file.endswith('.gz'):
            file = gzip.open(self.input_file, 'rt')
        else:
            file = open(self.input_file, 'r', encoding="UTF-8")

        while True:
            metadata_line = file.readline()
            sequence_lines = file.readline()
            if not metadata_line or not sequence_lines:
                break
            while True:
                pos = file.tell()
                line_next = file.readline()
                if not line_next:
                    break
                if line_next.strip() == "":
                    break
                if line_next.startswith('>'):
                    file.seek(pos)
                    break
                else:
                    sequence_lines += line_next.strip()

            yield metadata_line, sequence_lines

    def write_output(self, lock: Lock, read_metadata, forward_match, backward_match, inter_primer_sequence, possible_orf):
        with lock, open(self.output_file, 'a') as out_file:
            out_file.write(read_metadata.replace('|', ';').replace(',', ';').strip()
                           + f"{forward_match.score};{forward_match.read};{forward_match.start_index};"
                           + f"{backward_match.score};{backward_match.read};{backward_match.start_index};"
                           + f"{inter_primer_sequence};{str(possible_orf)}\n")
