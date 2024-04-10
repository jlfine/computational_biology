# Shuffle the codon content of any ORF
# Jacob L. Fine
# November 27th, 2023

import numpy as np
import pandas as pd

seq_file = pd.read_csv('C:\\Users\\jacob\\OneDrive\\U of T 2022-2023\\Blencowe\\Nov_2023_Blencowe\\counts_file_with_ORF_analysis_v5_first10lines.csv')

code = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'}


# human codon frequencies obtained from hg38 http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/Hsapi38-summary-codon.html
human_frqs_hg38 = {'GCT': 1.84,
'GCC': 2.77,
'GCG': 0.74,
'GCA': 1.58,
'GGT': 1.08,
'GGC': 2.22,
'GGG': 1.65,
'GGA': 1.65,
'CCT': 1.75,
'CCC': 1.98,
'CCG': 0.69,
'CCA': 1.69,
'ACT': 1.31,
'ACC': 1.89,
'ACG': 0.61,
'ACA': 1.51,
'GTT': 1.10,
'GTC': 1.45,
'GTG': 2.81,
'GTA': 0.71,
'TCT': 1.52,
'TCC': 1.77,
'TCG': 0.44,
'TCA': 1.22,
'AGT': 1.21,
'AGC': 1.95,
'CGT': 0.45,
'CGC': 1.04,
'CGG': 1.14,
'CGA': 0.62,
'AGG': 1.20,
'AGA': 1.22,
'CTT': 1.32,
'CTC': 1.96,
'CTG': 3.96,
'CTA': 0.72,
'TTG': 1.29,
'TTA': 0.77,
'TTT': 1.76,
'TTC': 2.03,
'AAT': 1.70,
'AAC': 1.91,
'AAG': 3.19,
'AAA': 2.44,
'GAT': 2.18,
'GAC': 2.51,
'GAG': 3.96,
'GAA': 2.90,
'CAT': 1.09,
'CAC': 1.51,
'CAG': 3.42,
'CAA': 1.23,
'ATT': 1.60,
'ATC': 2.08,
'ATA': 0.75,
'ATG': 2.20,
'TAT': 1.22,
'TAC': 1.53,
'TGA': 0.16,
'TAG': 0.08,
'TAA': 0.10,
'TGT': 1.06,
'TGC': 1.26,
'TGG': 1.32}


human_frqs_hg38 = {key: round(frq/100,4) for key, frq in human_frqs_hg38.items()}

print(human_frqs_hg38)
# human tRNA frequencies obtained from hg38 http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/Hsapi38-summary-codon.html
human_tRNA_copy_hg38 = {'AGC': 26,
'GGC': 0,
'CGC': 4,
'TGC': 8,
'ACC': 0,
'GCC': 14,
'CCC': 5,
'TCC': 9,
'AGG': 9,
'GGG': 0,
'CGG': 4,
'TGG': 7,
'AGT': 9,
'GGT': 0,
'CGT': 5,
'TGT': 6,
'AAC': 9,
'GAC': 0,
'CAC': 13,
'TAC': 5,
'AGA': 9,
'GGA': 0,
'CGA': 4,
'TGA': 4,
'ACT': 0,
'GCT': 8,
'ACG': 7,
'GCG': 0,
'CCG': 4,
'TCG': 6,
'CCT': 5,
'TCT': 6,
'AAG': 9,
'GAG': 0,
'CAG': 9,
'TAG': 3,
'CAA': 6,
'TAA': 4,
'AAA': 0,
'GAA': 10,
'ATT': 0,
'GTT': 25,
'CTT': 15,
'TTT': 12,
'ATC': 0,
'GTC': 13,
'CTC': 8,
'TTC': 8,
'ATG': 0,
'GTG': 9,
'CTG': 13,
'TTG': 6,
'AAT': 15,
'GAT': 3,
'TAT': 5,
'CAT': 10,
'ATA': 0,
'GTA': 13,
'CTA': 0,
'TTA': 0,
'ACA': 0,
'GCA': 29,
'CCA': 7,
'TCA': 1}


class GeneticCodeFreqs:
    def __init__(self, code, frqs):
        self.code = code
        self.frqs = human_frqs_hg38

    def get_amino_codon_freq(self):
        # initialize empty dictionary for posterior codon probabilities
        amino_codon_freq_dict = {}

        # go through the genetic code dictionary to initialize first layer of nested dict
        for codon, amino_acid in self.code.items():
            # adds each aa to the first layer of the dict
            if amino_acid not in amino_codon_freq_dict:
                amino_codon_freq_dict[amino_acid] = {}

            # looks up the codon frequency for the codon associated with that aa in the genetic code dict
            amino_codon_freq_dict[amino_acid][codon] = self.frqs[codon]

        # normalizes the frequencies for each codon given an aa by the freq of that amino acid
        for amino_acid, codons in amino_codon_freq_dict.items():
            total = sum(codons.values())  # total aa freq, associated with all the codons for that aa
            for codon, freq in codons.items():  # normalizes the codon frequency by the aa frequency
                amino_codon_freq_dict[amino_acid][codon] = round(freq / total,4)

        return amino_codon_freq_dict  # returns the amino_codon_freq_dict

# initializes the class with the genetic code and the set of human genomic codon frequencies
genetic_code = GeneticCodeFreqs(code, human_frqs_hg38)

# stores the nested posterior probability dict as a variable
amino_codon_freq_dict = genetic_code.get_amino_codon_freq()

# prints the nested dict
print(amino_codon_freq_dict)

seqs = seq_file['Sequence']
print(seqs)

