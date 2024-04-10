# Jacob L. Fine
# April 10th, 2024

# Finds the longest stretch of specific codons in an open reading frame (ORF).
# For instance, one may ask what the longest substring made of GGG, CCC, AAA codons is within a particular string.

codon_list = ['AAA','TTT','CCC','GGG']  # the list of codons we may desire to have in our sequence

# the minimal size of codon windows we choose to search, in codons
k = 1

# the percent of codons in the final stretch we wish to consist of our codon of interest, i.e.,
# p = 90% or more of our codons must be AAA, TTT, CCC, or GGG
p = 100


# the function we wish to use that finds the longest codon substring
def codon_stretch(seq, codon_list,k,p) -> str: 

    seq = [seq[i:i + 3] for i in range(0, len(seq), 3)]
    # print(seq)
    target_p = p  # this determines the amount of other codons (not in codon_list) that we allow into our stretch
    for i in range(len(seq), 0, - 1):
        for j in range(len(seq) - i + 1):  # goes through the codons
            substring = seq[j:j + i] # slices the string to obtain each substring
            if len(substring) < k:  # if the substring is too short, we set it to ''
                substring = ''
                coords = ''
                print(f'no substring for sequence')
                return substring, coords
            percentage = 0

            percentage = sum(codon in codon_list for codon in substring) / len(substring)

            if percentage >= target_p / 100:
                substring = "".join(substring) # joins the substirng which was originally a list
                coords = str(3*j) + '-' + str(len(substring)) # the coords of the match
                return substring, coords
            
seq = 'TATTATGGGCCCGGGCCCGGGCCCTATTATTATGGGCCCGGGCCC'  # the string we are searching

substring, coords = codon_stretch(seq, codon_list,k,p)
print(substring)  # prints the longest substring consistent of p percent of our codons of interest.
print(coords)  # the coords of the subtring in the string.


