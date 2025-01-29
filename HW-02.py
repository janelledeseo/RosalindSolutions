"""
Counting DNA Nucleotides
URL: https://rosalind.info/problems/dna/

Given: A DNA string s of length at most 1000 nt
Return: Four integers (separated by spaces) counting the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur in s
"""
def count_DNA(string):
    countA = string.count("A")
    countC = string.count("C")
    countG = string.count("G")
    countT = string.count("T")
    return countA, countC, countG, countT



"""
Transcribing DNA into RNA
URL: https://rosalind.info/problems/rna/

Given: A DNA string t having length at most 1000 nt
Return: The transcribed RNA string of t
"""
def dna2rna(string):
    return string.replace('T','U')



"""
Complementing a Strand of DNA
URL: https://rosalind.info/problems/revc/

Given: A DNA string s of length at most 1000 bp
Return: The reverse complement of sc of s
"""



"""
Computing GC Content 
URL: https://rosalind.info/problems/gc/

Given: At most 10 DNA strings in FASTA format (of length at most 1 kbp each)
Return: The ID of the string having the highest GC-content, followed by the GC-content of that string
"""



"""
Counting Point Mutations
URL: https://rosalind.info/problems/hamm/

Given: Two DNA strings s and t of equal length (not exceeding 1 kbp)
Return: The Hamming distance dH(s,t)
"""



"""
Translating RNA into Protein
URL: https://rosalind.info/problems/prot/

Given: An RNA string s corresponding to a strand of mRNA (of length at most 10 kbp)
Return: The protein string encoded by s
"""
