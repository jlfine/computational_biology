# Reverse translates every ORF in the list provided
# Jacob L. Fine
# January 11th, 2024


import random
import pandas as pd



seeds = [1001, 2002, 3003, 4004, 5005] # a list of seeds for different initializations, can toggle

s = seeds[2] # selects which seed, can toggle

random.seed(s)

seq_col = 'Sequence' # ensure the sequence column is 'Sequence' and contains each ORF

seq_file = pd.read_csv('/path/to/file.csv') # substitute with CSV of interest

# the genetic code dict
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

seqs = seq_file[seq_col]  # gets the seqs as a seq col
print(seqs)


class ORF_translator:
    def __init__(self,seq_file,code):
        self.seq_file = seq_file
        self.code = code

    def translate_orf(self, orf):
        polypeptide = ''
        for i in range(0, len(orf)-2, 3):
            codon = orf[i:i + 3]  # gets each codon in the ord
            if len(codon) == 3:
                if 'N' in codon:
                    polypeptide += 'X'  # adds X for codons that have 'N' in them

                else:
                    polypeptide += self.code[codon]  # looks up amino acid for each codon and appends it to polypeptide
            else:
                pass
        return polypeptide  # returns polypeptide for the row

    def add_polypeptide_column(self):
        self.seq_file['Polypeptide'] = self.seq_file['Sequence'].apply(self.translate_orf)  # applies the function to each row in the df



translator = ORF_translator(seq_file, code)
translator.add_polypeptide_column()


class ORF_reverse_translator:  # uses Monte Carlo simulation to randomly generate N possible ORFs for each peptide
    def __init__(self,seq_file,code):
        self.seq_file = seq_file
        self.code = code

    def reverse_translator(self,polypeptide):
        random_orf = ''
        for i in range(0, len(polypeptide), 1):
            amino_acid = polypeptide[i:i + 1]  # gets each codon in the ord
            if amino_acid != 'X':  # only considers defined amino acids
                possible_codon_dict_for_aa = amino_codon_freq_dict[amino_acid]  # gets possible codons dict for amino acid, and their probs
                possible_codon = random.choices(list(possible_codon_dict_for_aa.keys()),
                                                weights=list(possible_codon_dict_for_aa.values()),k=1)  # randomly selects a codon based on the posterior probabilities, given the amino acid, from the dict
                random_orf += possible_codon[0]
        return random_orf  # returns a random possible for the row

    def add_random_orf(self):
        self.seq_file['Sequence_Random'] = self.seq_file['Polypeptide'].apply(self.reverse_translator)  # applies the function to each row in the df



reverse_translator = ORF_reverse_translator(seq_file, code)
reverse_translator.add_random_orf()

seq_file.to_csv(f'/path/to/file_output_seed={s}.csv')
