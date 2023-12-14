"""
DNA Processing Script with Argument Parsing

Benjamin Rudski

Intermediate Python Skills Workshop, MiCM - Fall 2023
"""

import argparse
import dataclasses
import os
import shutil
from typing import Dict, List

import numpy as np

class RnaSequence:
    """
    RNA Sequence

    This class represents an mRNA sequence.

    Attributes:
        - rna_sequence: string of ribonucleotides (A, U, C, G).
        - codon_table: dictionary containing codons as keys and amino acids as values.

    Methods:
        - translate: convert mRNA into amino acid sequence.
    """

    rna_sequence: str
    codon_table: Dict[str, str]

    def __init__(self, rna_sequence: str):
        self.rna_sequence = rna_sequence

        amino_acid_to_codon_table = {
            "F": ["UUU", "UUC"],
            "L": ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
            "I": ["AUU", "AUC", "AUA"],
            "M": ["AUG"],
            "V": ["GUU", "GUC", "GUA", "GUG"],
            "S": ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],
            "P": ["CCU", "CCC", "CCA", "CCG"],
            "T": ["ACU", "ACC", "ACA", "ACG"],
            "A": ["GCU", "GCC", "GCA", "GCG"],
            "Y": ["UAU", "UAC"],
            "STOP": ["UAA", "UAG", "UGA"],
            "H": ["CAU", "CAC"],
            "Q": ["CAA", "CAG"],
            "N": ["AAU", "AAC"],
            "K": ["AAA", "AAG"],
            "D": ["GAU", "GAC"],
            "E": ["GAA", "GAG"],
            "C": ["UGU", "UGC"],
            "W": ["UGG"],
            "R": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
            "G": ["GGU", "GGC", "GGA", "GGG"]
        }

        # Convert the table to have the codons as keys and the amino acids as values.
        self.codon_table = {}

        for amino_acid, codon_list in amino_acid_to_codon_table.items():
            for codon in codon_list:
                self.codon_table[codon] = amino_acid

    def __str__(self) -> str:
        return f"RNA sequence of {len(self.rna_sequence)} nucleotides"

    def translate(self) -> str:
        """
        Translate mRNA to Amino Acid Sequence.

        Translate mRNA sequence containing a start codon and open reading frame
        to an amino acid sequence using single-letter codes.
        
        Returns:
            - string of single-letter amino acids, beginning with M.
        """

        amino_acid_sequence = ""
        start_position = -1

        # We want to make sure that we don't go beyond the end of the string!
        for i in range(len(self.rna_sequence)-2):

            candidate_codon = self.rna_sequence[i: i+3]

            if candidate_codon == "AUG":
                start_position = i
                amino_acid_sequence += self.codon_table["AUG"]
                break

        # Again, we want to not go beyond the end.
        for i in range(start_position + 3, len(self.rna_sequence)-2, 3):
            
            new_codon = self.rna_sequence[i: i+3]
            new_amino_acid = self.codon_table[new_codon]

            if new_amino_acid == "STOP":
                break
            else:
                amino_acid_sequence += new_amino_acid

        return amino_acid_sequence


@dataclasses.dataclass
class DnaSequence:
    """
    DNA Sequence

    Data class representing a DNA sequence.

    Attributes:
        - dna_sequence: string of nucleotides.
        - sequence_name: name from FASTA file.

    Methods:
        - transcribe: transcribe DNA to mRNA
    """

    dna_sequence: str
    sequence_name: str

    def transcribe(self, is_template_strand: bool = True) -> RnaSequence:
        """
        Transcribe DNA to mRNA.

        ...
        """
        
        if is_template_strand:
            m_rna_sequence = ""

            for nt in self.dna_sequence:
                if nt == "A":
                    m_rna_sequence += "U"
                elif nt == "T":
                    m_rna_sequence += "A"
                elif nt == "C":
                    m_rna_sequence += "G"
                else:
                    m_rna_sequence += "C"
        else:
            m_rna_sequence = self.dna_sequence.replace("T", "U")

        return RnaSequence(m_rna_sequence)


def compute_cg_proportion(sequences: List[DnaSequence]) -> np.ndarray:
    """
    Compute the CG percentage at each position in sequences of the same length.

    This function computes the percentage of nucleotides that are C or G at
    each position in a set of nucleotide sequences of the same length.

    Parameters:
        - sequences: List containing DNA sequences of the **same length**.

    Returns:
        - NumPy array of shape (sequence_length,) where sequence_length 
          corresponds to the common length of the nucleotide sequences.
    """

    # Split sequences into the nucleotides (see https://www.geeksforgeeks.org/python-split-string-into-list-of-characters/)
    split_sequences = [list(seq.dna_sequence.upper()) for seq in sequences]

    # Create NumPy array based on the sequences
    nucleotide_array = np.array(split_sequences)

    # Create the voting array, with the same shape as the nucleotide array
    cg_count_array = np.zeros(nucleotide_array.shape)

    # Now, for the assignment magic! We're using conditional assignment!
    cg_count_array[nucleotide_array=="C"] = 1
    cg_count_array[nucleotide_array=="G"] = 1

    # Now, for the mean
    cg_probs = cg_count_array.mean(axis=0)

    return cg_probs


def read_dna_fasta_file(filename: str, max_sequences: int = 0) -> List[DnaSequence]:
    """
    Read DNA FASTA file.

    This function reads a DNA FASTA file and creates a list of ``DnaSequence`` objects.

    Parameters:
        - filename: path to the FASTA file.
        - max_sequences: Maximum number of sequences to process in the file. If zero is passed, all sequences will be read.

    Returns:
        - List of ``DnaSequence`` objects representing the loaded DNA sequences.
    """

    my_sequences: List[DnaSequence] = []

    with open(filename) as fasta_file:

        current_sequence = ""
        current_name = ""

        for line in fasta_file:

            current_line = line.strip()

            if current_line[0] == ">":

                if current_sequence != "":
                    new_sequence_object = DnaSequence(current_sequence, current_name)
                    my_sequences.append(new_sequence_object)

                    current_sequence = ""

                sequence_name = current_line[1:].strip()
                current_name = sequence_name

                # End if we've hit the maximum number of sequences
                if max_sequences > 0 and len(my_sequences) >= max_sequences:
                    break

            else:
                current_sequence += current_line

    return my_sequences


def write_fasta_file(sequences: List[str], names: List[str], filename: str):
    """
    Write FASTA file containing sequences of any type.

    This function writes a FASTA file, regardless of the type of sequences passed in.
    Lines in the FASTA file have a fixed width of 80 characters.

    Parameters:
        - sequences: sequences to write. May be nucleotide or amino acids.
        - names: sequence names, or other metadata to encode in the FASTA file.
        - filename: path to the FASTA file.

    Returns:
        - None.
    """

    with open(filename, "w") as fasta_file:
        for seq, name in zip(sequences, names):
            fasta_file.write(f"> {name}\n")

            # We want to clip the lines to be a certain length so that the file is easier to read.
            # Let's say 80 characters.
            line_length = 80
            seq_length = len(seq)

            number_of_full_lines = seq_length // line_length
            remaining_chars = seq_length % line_length

            for i in range(number_of_full_lines):
                
                # Get the sub_sequence
                sub_sequence = seq[i * line_length: (i+1) * line_length]
                fasta_file.write(f"{sub_sequence}\n")

            # Check if we have a remainder
            if remaining_chars > 0:
                # Take the last `remaining_chars` characters
                sub_sequence = seq[-remaining_chars:]
                fasta_file.write(f"{sub_sequence}\n")


def process_fasta_file(fasta_file: str, output_dir: str, is_template_strand: bool = True, compute_stats: bool = False, max_sequences: int = 0):
    """
    FASTA File Processing

    Translate DNA from a specified FASTA file to a new FASTA file containing the protein sequences.

    Parameters:
        - fasta_file: FASTA file containing the DNA sequences.
        - output_dir: Name for output directory where the proteins are stored.
        - is_template_strand: Indicate whether DNA is on template strand.
        - compute_stats: Indicate whether the GC stats should be computed and saved.
        - max_sequences: Maximum number of sequences to process in the file. If zero is passed, all sequences will be analysed.
    """

    # Create our output directory, if necessary.
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir) 

    # Perform the translation
    dna_sequences = read_dna_fasta_file(fasta_file, max_sequences=max_sequences)
    protein_sequences: List[str] = []
    protein_sequence_names: List[str] = []

    for sequence in dna_sequences:
        m_rna_sequence = sequence.transcribe(is_template_strand=is_template_strand)
        protein_sequence = m_rna_sequence.translate()

        protein_sequences.append(protein_sequence)
        protein_sequence_names.append(sequence.sequence_name)

    original_dna_base_filename = os.path.basename(fasta_file)
    original_dna_base_name_no_ext, original_dna_base_name_ext = os.path.splitext(original_dna_base_filename)

    protein_filename = f"{original_dna_base_name_no_ext}_PROTEIN{original_dna_base_name_ext}"

    final_protein_filename = os.path.join(output_dir, protein_filename)
    
    # Write the FASTA file
    write_fasta_file(sequences=protein_sequences, names=protein_sequence_names, filename=final_protein_filename)

    # Compute the stats, if requested
    if compute_stats:
        stats_path = os.path.join(output_dir, "Stats")

        if not os.path.isdir(stats_path):
            os.mkdir(stats_path)
   

        stats_array = compute_cg_proportion(sequences=dna_sequences)

        stats_filename = f"{original_dna_base_name_no_ext}_STATS.txt"
        full_stats_filename = os.path.join(stats_path, stats_filename)

        np.savetxt(full_stats_filename, stats_array)


def run_overall_pipeline(sequence_folder: str, output_folder: str, compute_stats: bool = True, max_sequences: int = 0):
    """
    Run the overall pipeline on a set of FASTA files.

    This function takes all FASTA files from a specified folder, performs transcription
    and translation, and computes the position-wise CG proportions for all sequences in
    a given file.

    Parameters:
        - sequence_folder: Folder containing FASTA files to process.
        - output_folder: Folder where the output proteins and statistics will be stored.
        - compute_stats: Indicate whether the CG proportions should be computed.
        - max_sequences: Maximum number of sequences to process in each file. If zero is passed, all sequences will be analysed.
    """

    # Iterate through all files in the sequence folder.

    files_to_consider: List[str] = [f for f in os.listdir(sequence_folder) if f.endswith(".fasta")]

    for fasta_file in files_to_consider:
        is_template_strand = not "non_template" in fasta_file
        full_input_path = os.path.join(sequence_folder, fasta_file)
        process_fasta_file(fasta_file=full_input_path, output_dir=output_folder, compute_stats=compute_stats, is_template_strand=is_template_strand, max_sequences=max_sequences)


if __name__ == "__main__":
    
    # Create the parser
    parser = argparse.ArgumentParser(
        description="Batch DNA Processor",
        epilog="MiCM Intermediate Python Workshop -- Fall 2023"
    )

    # Add the arguments
    parser.add_argument("input_dir", type=str, help="Directory containing FASTA files.")
    parser.add_argument("-o", "--output-dir", type=str, default="./protein_output", help="Output directory for protein FASTA files. Default: `protein_output` in the current directory.")
    parser.add_argument("-m", "--max-sequences", type=int, default=0, help="Maximum number of sequences to consider (0=all).")
    parser.add_argument("--stats", action="store_true", help="Compute and save the GC fraction statistics.")

    # Parse the arguments
    args = parser.parse_args()

    # Retrieve the arguments using the names we passed in above when adding them.
    sequence_folder = args.input_dir
    output_folder = args.output_dir
    max_sequences = args.max_sequences
    compute_stats = args.stats

    run_overall_pipeline(sequence_folder=sequence_folder, output_folder=output_folder, compute_stats=compute_stats, max_sequences=max_sequences)

