# SUBS/PRTM/SPLC/REVP/TRAN/LCSM/ORF/

"""
Finding a Motif in DNA
URL: http://rosalind.info/problems/subs/

Given: Two DNA strings s and t (each of length at most 1 kbp).
Return: All locations of t as a substring of s.
"""
def subs(string1, string2):
    loc = []
    for i in range(len(string1)):
        if string2 == string1[i: i+len(string2)]:
            loc.append(i+1)
    return loc



"""
Calculating Protein Mass
URL: http://rosalind.info/problems/prtm/

Given: A protein string P of length at most 1000 aa.
Return: The total weight of P. Consult the monoisotopic mass table.
"""
mass = {
    'A': 71.03711,  'C': 103.00919, 'D': 115.02694,
    'E': 129.04259, 'F': 147.06841, 'G': 57.02146,
    'H': 137.05891, 'I': 113.08406, 'K': 128.09496,
    'L': 113.08406, 'M': 131.04049, 'N': 114.04293,
    'P': 97.05276,  'Q': 128.05858, 'R': 156.10111,
    'S': 87.03203,  'T': 101.04768, 'V': 99.06841,
    'W': 186.07931, 'Y': 163.06333
}

def protein_mass(P_string):
    return sum(mass[i] for i in P_string if i in mass)



"""
RNA Splicing
URL: http://rosalind.info/problems/splc/

Given: A DNA string s (of length at most 1 kbp) and a collection of substrings of s acting as introns. All strings are given in FASTA format.
Return: A protein string resulting from transcribing and translating the exons of s. (Note: Only one solution will exist for the dataset provided.)
"""
coding_table = {
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 
    'UUU': 'F', 'UUC': 'F',
    'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'UAU': 'Y', 'UAC': 'Y',
    'UGU': 'C', 'UGC': 'C',
    'UGG': 'W', 
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H',
    'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I',
    'AUG': 'M', 
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 
    'AAU': 'N', 'AAC': 'N', 
    'AAA': 'K', 'AAG': 'K',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAU': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    'UAA': '*', 'UAG': '*', 'UGA': '*'  # Stop codons
}

def RNA_splicing(dna_string, introns):
    # Remove introns
    for intron in introns:
        dna_string = dna_string.replace(intron, "")

    # Transcribe DNA to RNA
    rna_string = dna_string.replace("T", "U")

    # Translate RNA to Protein
    protein_string = ""
    for i in range(0, len(rna_string) - 2, 3):
        codon = rna_string[i:i+3]
        amino_acid = coding_table.get(codon, '')
        protein_string += amino_acid
        if amino_acid == '*':
            break

    return protein_string



"""
Locating Restriction Sites
URL: http://rosalind.info/problems/revp/

Given: A DNA string of length at most 1 kbp in FASTA format.
Return: The position and length of every reverse palindrome in the string having length between 4 and 12. You may return these pairs in any order.
"""
# Return the reverse complement of a DNA sequence
def switch(string):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(string))

# Find and print reverse palindromic sequences of length 4 to 12
def palindrome(string):
    for i in range(len(string)):
        for j in range(4, 13):  
            if i + j <= len(string):  # Ensure valid substring
                substring = string[i:i+j]
                if substring == switch(substring):  # Check if reverse complement matches
                    print(i + 1, j)  # Output 1-based position and length



"""
Transitions and Transversions
URL: http://rosalind.info/problems/tran/

Given: Two DNA strings s1 and s2 of equal length (at most 1 kbp).
Return: The transition/transversion ratio R(s1,s2).
"""
transitions = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
transversions = {('A', 'T'), ('T', 'A'), ('A', 'C'), ('C', 'A'),
                 ('T', 'G'), ('G', 'T'), ('C', 'G'), ('G', 'C')}

transition_c, transversions_c = 0, 0 # counters

for a, b in zip(s1, s2):
    if (a, b) in transitions:
        transition_c += 1
    elif (a, b) in transversions:
        transversions_c += 1

if transversions_c == 0:
    print("No transversions found, ratio is undefined")
else:
    print(transition_c / transversions_c)



"""
Finding a Shared Motif
URL: http://rosalind.info/problems/lcsm/

Given: A collection of k (k≤100) DNA strings of length at most 1 kbp each in FASTA format.
Return: A longest common substring of the collection. (If multiple solutions exist, you may return any single solution.)
"""



"""
Open Reading Frames
URL:http://rosalind.info/problems/orf/

Given: A DNA string s of length at most 1 kbp in FASTA format.
Return: Every distinct candidate protein string that can be translated from ORFs of s. Strings can be returned in any order.
"""
