# Jacob L. Fine
# Aug 23, 2023
# Description: analyze ORFs of the genes of interest


import pandas as pd
import numpy as np
from statistics import geometric_mean
import math
from collections import Counter
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--dir", type=str, help='the directory of files')
parser.add_argument("--file", type=str, help="the seq file")
parser.add_argument("--col_name", type=str, help="the seq col")
parser.add_argument("--short_name", type=str, help="the brief name of file")

args = parser.parse_args()

directory = str(args.dir)
file_name = str(args.file)
col_name = str(args.col_name)
short_name = str(args.short_name)

print(col_name)

merged_df_with_seqs = pd.read_csv(file_name)  # the file with all the seqs



# tAI weights derived from weighted sum of tRNA copy number according to https://www.ncbi.nlm.nih.gov/pmc/articles/PMC521650/
codon_W_dict_human = {'ATA': 0.588322, 'ATC': 0.213397, 'ATT': 0.481031, 'ATG': 0.296384, 'ACA': 0.441909,
                      'ACC': 0.074689, 'ACG': 0.269117, 'ACT': 0.266746, 'AAC': 0.74096, 'AAT': 0.303794,
                      'AAA': 0.355661, 'AAG': 0.686426, 'AGC': 0.237107, 'AGT': 0.097214, 'AGA': 0.17783,
                      'AGG': 0.269117, 'CTA': 0.352993, 'CTC': 0.074689, 'CTG': 0.327208, 'CTT': 0.266746,
                      'CCA': 0.471547, 'CCC': 0.074689, 'CCG': 0.259632, 'CCT': 0.266746, 'CAC': 0.266746,
                      'CAT': 0.109366, 'CAA': 0.17783, 'CAG': 0.506224, 'CGA': 0.383225, 'CGC': 0.058091,
                      'CGG': 0.239478, 'CGT': 0.207469, 'GTA': 0.41227, 'GTC': 0.074689, 'GTG': 0.48607,
                      'GTT': 0.266746, 'GCA': 1.0, 'GCC': 0.215768, 'GCG': 0.279787, 'GCT': 0.770599,
                      'GAC': 0.385299, 'GAT': 0.157973, 'GAA': 0.237107, 'GAG': 0.39834, 'GGA': 0.266746,
                      'GGC': 0.414938, 'GGG': 0.329579, 'GGT': 0.170124, 'TCA': 0.382632, 'TCC': 0.074689,
                      'TCG': 0.19917, 'TCT': 0.266746, 'TTC': 0.296384, 'TTT': 0.121517, 'TTA': 0.118554,
                      'TTG': 0.258447, 'TAC': 0.385299, 'TAT': 0.157973, 'TGC': 0.859514, 'TGT': 0.352401,
                      'TGG': 0.227623}

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
    'TAC': 'Y', 'TAT': 'Y', 'TAA': ')', 'TAG': ')',
    'TGC': 'C', 'TGT': 'C', 'TGA': ')', 'TGG': 'W', 'NNN': 'X'}  # accounts for unknown codons


def tAI_human_f0(sequence):
    if len(sequence) >= 3 and ('Seq' not in sequence):  # bigger than or equal to 3 st it will have 1 codon at least
        seq_con = [sequence[i:i + 3] for i in range(0, len(sequence), 3) if len(sequence[
                                                                                i:i + 3]) % 3 == 0]  # converst seq to list   # converts the seq into a list so it only takes the frame
        # print(seq_con)
        if not (len(seq_con) == 1 and seq_con[0] in ['TGA', 'TAA', 'TAG']):  # if there is not just one stop codon
            tAI_seq = geometric_mean([codon_W_dict_human[codon] for codon in seq_con if (
                    (codon not in ['TGA', 'TAA', 'TAG']) and ('N' not in codon) and (len(codon) % 3 == 0))])
            return tAI_seq
    else:
        return np.nan


def tAI_human_f1(sequence):
    sequence = sequence[1:]  # does a frame shift by 1
    if len(sequence) >= 3 and ('Seq' not in sequence):  # bigger than or equal to 3 st it will have 1 codon at least
        seq_con = [sequence[i:i + 3] for i in range(0, len(sequence), 3) if len(sequence[
                                                                                i:i + 3]) % 3 == 0]  # converst seq to list   # converts the seq into a list so it only takes the frame
        # print(seq_con)
        if not (len(seq_con) == 1 and seq_con[0] in ['TGA', 'TAA', 'TAG']):  # if there is not just one stop codon
            tAI_seq = geometric_mean([codon_W_dict_human[codon] for codon in seq_con if (
                    (codon not in ['TGA', 'TAA', 'TAG']) and ('N' not in codon) and (len(codon) % 3 == 0))])
            return tAI_seq
    else:
        return np.nan


def tAI_human_f2(sequence):
    sequence = sequence[2:]  # does a frame shift by 1
    if len(sequence) >= 3 and ('Seq' not in sequence):  # bigger than or equal to 3 st it will have 1 codon at least
        seq_con = [sequence[i:i + 3] for i in range(0, len(sequence), 3) if len(sequence[
                                                                                i:i + 3]) % 3 == 0]  # converst seq to list   # converts the seq into a list so it only takes the frame
        # print(seq_con)
        if not (len(seq_con) == 1 and seq_con[0] in ['TGA', 'TAA', 'TAG']):  # if there is not just one stop codon
            tAI_seq = geometric_mean([codon_W_dict_human[codon] for codon in seq_con if (
                    (codon not in ['TGA', 'TAA', 'TAG']) and ('N' not in codon) and (len(codon) % 3 == 0))])
            return tAI_seq
    else:
        return np.nan


def shannon_entropy_nuc(sequence):
    if len(sequence) > 0 and ('Seq' not in sequence):
        # Calculate the frequency of each nucleotide in the sequence
        nucleotide_counts = Counter(sequence)
        # Calculate the probability of each nucleotide
        nucleotide_probabilities = [count / len(sequence) for count in nucleotide_counts.values()]
        # Calculate the Shannon entropy
        # print(nucleotide_probabilities)
        entropy = -sum(p * math.log2(p) for p in nucleotide_probabilities)
        entropy = round(entropy, 4)

        return entropy
    else:
        return np.nan


def shannon_entropy_codon(sequence):
    if len(sequence) >= 3 and ('Seq' not in sequence):  # bigger than or equal to 3 st it will have 1 codon at least
        seq_con = [sequence[i:i + 3] for i in range(0, len(sequence), 3) if len(sequence[
                                                                                i:i + 3]) % 3 == 0]  # converst seq to list   # converts the seq into a list so it only takes the frame
        # print(seq_con)

        # Calculate the frequency of each nucleotide in the sequence
        codon_counts = Counter(seq_con)
        # Calculate the probability of each nucleotide
        codon_probabilities = [count / len(seq_con) for count in codon_counts.values()]
        # Calculate the Shannon entropy
        # print(nucleotide_probabilities)
        entropy = -sum(p * math.log2(p) for p in codon_probabilities)
        entropy = round(entropy, 4)

        return entropy
    else:
        return np.nan


def shannon_entropy_protein(sequence):
    if len(sequence) >= 3 and ('Seq' not in sequence):  # bigger than or equal to 3 st it will have 1 codon at least
        seq_con = [sequence[i:i + 3] for i in range(0, len(sequence), 3) if len(sequence[
                                                                                i:i + 3]) % 3 == 0]  # converst seq to list   # converts the seq into a list so it only takes the frame
        # print(seq_con)

        seq_con = ['NNN' if 'N' in codon else codon for codon in
                   seq_con]  # converges N containing codons to NNN so they can be translated to X

        seq_con = ([code[codon] for codon in seq_con if (
                (codon not in ['TGA', 'TAA', 'TAG']) and (len(codon) % 3 == 0))])
        # Calculate the frequency of each nucleotide in the sequence
        aa_counts = Counter(seq_con)
        # Calculate the probability of each nucleotide
        aa_probabilities = [count / len(seq_con) for count in aa_counts.values()]
        # Calculate the Shannon entropy
        # print(nucleotide_probabilities)
        entropy = -sum(p * math.log2(p) for p in aa_probabilities)
        entropy = round(entropy, 4)
        return entropy
    else:
        return np.nan


def GC(sequence):
    if len(sequence) > 0 and ('Seq' not in sequence):
        # Count the number of G and C nucleotides in the sequence
        gc_count = sequence.count('G') + sequence.count('C')
        # Calculate the GC content as a fraction of the total length of the sequence
        gc_fraction = round(gc_count / len(sequence), 3)
        return gc_fraction
    else:
        return np.nan


def AG(sequence):
    if len(sequence) > 0 and ('Seq' not in sequence):
        # Count the number of G and C nucleotides in the sequence
        ag_count = sequence.count('A') + sequence.count('G')
        # Calculate the GC content as a fraction of the total length of the sequence
        ag_fraction = round(ag_count / len(sequence), 3)
        return ag_fraction
    else:
        return np.nan


def CU(sequence):
    if len(sequence) > 0 and ('Seq' not in sequence):
        # Count the number of G and C nucleotides in the sequence
        ag_count = sequence.count('C') + sequence.count('T')
        # Calculate the GC content as a fraction of the total length of the sequence
        ag_fraction = round(ag_count / len(sequence), 3)
        return ag_fraction
    else:
        return np.nan


def A(sequence):
    if len(sequence) > 0 and ('Seq' not in sequence):
        # Count the number of G and C nucleotides in the sequence
        ag_count = sequence.count('A')
        # Calculate the GC content as a fraction of the total length of the sequence
        ag_fraction = round(ag_count / len(sequence), 3)
        return ag_fraction
    else:
        return np.nan


def G(sequence):
    if len(sequence) > 0 and ('Seq' not in sequence):
        # Count the number of G and C nucleotides in the sequence
        ag_count = sequence.count('G')
        # Calculate the GC content as a fraction of the total length of the sequence
        ag_fraction = round(ag_count / len(sequence), 3)
        return ag_fraction
    else:
        return np.nan


def U(sequence):
    if len(sequence) > 0 and ('Seq' not in sequence):
        # Count the number of G and C nucleotides in the sequence
        ag_count = sequence.count('T')
        # Calculate the GC content as a fraction of the total length of the sequence
        ag_fraction = round(ag_count / len(sequence), 3)
        return ag_fraction
    else:
        return np.nan


def C(sequence):
    if len(sequence) > 0 and ('Seq' not in sequence):
        # Count the number of G and C nucleotides in the sequence
        ag_count = sequence.count('C')
        # Calculate the GC content as a fraction of the total length of the sequence
        ag_fraction = round(ag_count / len(sequence), 3)
        return ag_fraction
    else:
        return np.nan


def GC3(sequence):
    if len(sequence) > 0 and ('Seq' not in sequence):
        seq_con = [sequence[i:i + 3] for i in range(0, len(sequence), 3) if len(sequence[
                                                                                i:i + 3]) % 3 == 0]  # converst seq to list   # converts the seq into a list so it only takes the frame
        seq_third_position = [codon[2] for codon in seq_con]  # only takes the third position of each codon
        seq_third_position_str = ''.join(seq_third_position)
        if len(seq_third_position_str) > 0:
            GC3_seq = round(
                (seq_third_position_str.count('G') + seq_third_position_str.count('C')) / len(seq_third_position_str),
                3)
        else:
            GC3_seq = np.nan
        return GC3_seq
    else:
        return np.nan


def ACU3(sequence):
    if len(sequence) > 0 and ('Seq' not in sequence):
        seq_con = [sequence[i:i + 3] for i in range(0, len(sequence), 3) if len(sequence[
                                                                                i:i + 3]) % 3 == 0]  # converst seq to list   # converts the seq into a list so it only takes the frame
        seq_third_position = [codon[2] for codon in seq_con]  # only takes the third position of each codon
        seq_third_position_str = ''.join(seq_third_position)
        if len(seq_third_position_str) > 0:
            ACU3_seq = round((seq_third_position_str.count('A') + seq_third_position_str.count(
                'T') + seq_third_position_str.count('C')) / len(seq_third_position_str), 3)
        else:
            ACU3_seq = np.nan
        return ACU3_seq
    else:
        return np.nan

def orf_length(sequence):
    length = len(sequence)
    return length

def rare_codon_count_0(sequence) -> int:  # maps str to int
    rare_codons = ['TCG', 'CGT', 'ACG', 'CGA', 'CCG', 'GTA', 'CTA', 'GCG', 'ATA',
                   'TTA']  # codons with less than 1% frequency in Hsapi38

    sequence = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]
    rare_count = sum(codon in rare_codons for codon in sequence)
    print(rare_count)
    return rare_count


def rare_codon_count_1(sequence) -> int:  # maps str to int
    sequence = sequence[1:]  # does a frame shift by 1
    rare_codons = ['TCG', 'CGT', 'ACG', 'CGA', 'CCG', 'GTA', 'CTA', 'GCG', 'ATA',
                   'TTA']  # codons with less than 1% frequency in Hsapi38

    sequence = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]
    rare_count = sum(codon in rare_codons for codon in sequence)
    return rare_count


def rare_codon_count_2(sequence) -> int:  # maps str to int
    sequence = sequence[2:]  # does a frame shift by 1
    rare_codons = ['TCG', 'CGT', 'ACG', 'CGA', 'CCG', 'GTA', 'CTA', 'GCG', 'ATA',
                   'TTA']  # codons with less than 1% frequency in Hsapi38

    sequence = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]
    rare_count = sum(codon in rare_codons for codon in sequence)
    return rare_count


def rare_codon_expected(sequence) -> int:  # maps str to int

    sequence = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]
    n_codons = len(sequence)
    p = 0.065  # percent of codons that are rare, with less than 1% frequency according to http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/Hsapi38-summary-codon.html
    n_expected = round(p * n_codons)
    return n_expected



# a list of all the functions to apply
function_list = [tAI_human_f0, tAI_human_f1, tAI_human_f2, rare_codon_count_0, rare_codon_count_1, rare_codon_count_2,
                 rare_codon_expected, orf_length, shannon_entropy_nuc, shannon_entropy_codon, shannon_entropy_protein,
                 A, G, C, U, AG, CU, GC, GC3, ACU3]
for function in function_list:
    merged_df_with_seqs[f'{function.__name__}'] = merged_df_with_seqs[col_name].apply(
        function)  # applies the functions to the column in a for loop  # gets the function name

merged_df_with_seqs.to_csv(directory + f'{short_name}.csv')
