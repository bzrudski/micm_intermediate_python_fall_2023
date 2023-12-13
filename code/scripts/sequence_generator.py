"""
Sequence Generator
Benjamin Rudski, 2023

This module contains functions to generate DNA sequences.

"""

import math
import random
import textwrap

NUCLEOTIDES = ['A', 'T', 'C', 'G']

NON_TEMPLATE_START_CODON = 'ATG'
NON_TEMPLATE_STOP_CODONS = ["TAA", "TAG", "TGA"]

TEMPLATE_START_CODON = 'TAC'
TEMPLATE_STOP_CODONS = ["ATT", "ATC", "ACT"]

def generate_dna_sequence(
        length: int, farthest_start: float = 0.1, 
        nt_weights: list[float] = [1, 1, 1, 1], 
        is_template_strand: bool = False
    ) -> str:
    """
    Generate random DNA sequence containing a gene.

    This function creates random sequences of DNA, containing a gene. This function is guaranteed to return
    an open reading frame that begins with AUG and ends with a stop codon.

    :param length: Length of the DNA sequence.
    :param farthest_start: farthest relative starting position within the sequence.
    :param nt_weights: weights of ``A``, ``T``, ``C`` and ``G`` nucleotides when constructing the sequence.
    :param is_template_strand: indicate whether the produced strand is the template strand.

    :return: the generated DNA sequence.
    """

    farthest_start_position = math.floor(length * farthest_start)

    start_position = random.randint(0, farthest_start_position)

    my_sequence = ""

    added_nucleotides = "".join(random.choices(population=NUCLEOTIDES, weights=nt_weights, k=start_position - 1))

    my_sequence += added_nucleotides

    if is_template_strand:
        my_sequence += TEMPLATE_START_CODON
    else:
        my_sequence += NON_TEMPLATE_START_CODON

    current_position = len(my_sequence) - 1

    farthest_stop_position = length - 3

    stop_position = random.randint(current_position + 3, farthest_stop_position)

    number_of_nucleotides_to_add = stop_position - current_position

    coding_nucleotides = "".join(random.choices(population=NUCLEOTIDES, weights=nt_weights, k=number_of_nucleotides_to_add))

    my_sequence += coding_nucleotides

    number_of_remaining_nucleotides_to_add = length - len(my_sequence)

    remaining_nucleotides = "".join(random.choices(population=NUCLEOTIDES, weights=nt_weights, k=number_of_remaining_nucleotides_to_add))

    my_sequence += remaining_nucleotides

    return my_sequence


def generate_random_sequences(
        output_file: str, header_prefix: str, 
        length: int = 1000, number_of_sequences: int = 100,
        nt_weights: list[float] = [1, 1, 1, 1],
        farthest_start: float = 0.2,
        fasta_line_length = 80,
        is_template_strand: bool = False):
    """
    Generate a series of random DNA sequences and save as FASTA.

    Generate a series of random coding DNA sequences of common length and save
    as a FASTA file. The metadata for the FASTA is specified in the ``header_prefix``
    argument. The number of the sequence is appended onto the header. Finally, text
    is added indicating whether the strand is template or non-template.

    :param output_file: output file for all the DNA sequences.
    :param header_prefix: label for the non-sequence lines (follows the ``>`` sign).
    :param length: length of each sequence.
    :param number_of_sequences: number of sequences to generate.
    :param fasta_line_length: length of the sequence lines in the FASTA file.
    :param farthest_start: farthest relative starting position within the sequence.
    :param nt_weights: weights of ``A``, ``T``, ``C`` and ``G`` nucleotides when constructing the sequence.
    :param is_template_strand: indicate whether the generated sequences are on the template strand.
    :return: None.
    """

    my_sequences = [
        generate_dna_sequence(
            length=length, 
            farthest_start=farthest_start,
            is_template_strand=is_template_strand,
            nt_weights=nt_weights
        )
        for _ in range(number_of_sequences)
    ]

    sequence_header_suffix = "TEMPLATE" if is_template_strand else "NON_TEMPLATE"
    sequence_header_template = f"> {header_prefix} {{}} {sequence_header_suffix}"

    my_sequence_headers = [
        sequence_header_template.format(i) for i in range(number_of_sequences)
    ]

    with open(output_file, "w") as fasta_file:
        for seq, header in zip(my_sequences, my_sequence_headers):
            fasta_file.write(header)
            fasta_file.write("\n")
            wrapped_sequence_lines = textwrap.wrap(seq, width=fasta_line_length)

            corrected_lines = [f"{line}\n" for line in wrapped_sequence_lines]

            fasta_file.writelines(corrected_lines)

            fasta_file.flush()

