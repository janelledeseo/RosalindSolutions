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
def rev_comp(string):
  complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
  reverse_complement = ''.join([complement[base] for base in string[::-1]])
  return reverse_complement



"""
Computing GC Content 
URL: https://rosalind.info/problems/gc/

Given: At most 10 DNA strings in FASTA format (of length at most 1 kbp each)
Return: The ID of the string having the highest GC-content, followed by the GC-content of that string
"""
def gc(strings):
    gc_contents = {}
    for k, v in strings.items():
        gc_content = (v.count("G") + v.count("C")) / len(v)
        gc_contents[k] = gc_content
    gc_contents = sorted(gc_contents.items(), key=lambda d: d[1], reverse=True) # key=lambda d: d[1] sorts dictionary items by their second element (ID, GC content)
    highest_gc_id, highest_gc_content = gc_contents[0]
    highest_gc_percentage = highest_gc_content * 100 
    return highest_gc_id, highest_gc_percentage



"""
Counting Point Mutations
URL: https://rosalind.info/problems/hamm/

Given: Two DNA strings s and t of equal length (not exceeding 1 kbp)
Return: The Hamming distance dH(s,t)
"""
def hamm(s, t):
    distance = 0
    assert len(s) == len(t)
    for i in range(len(s)):
        if s[i] != t[i]:
            distance += 1
    return distance



"""
Translating RNA into Protein
URL: https://rosalind.info/problems/prot/

Given: An RNA string s corresponding to a strand of mRNA (of length at most 10 kbp)
Return: The protein string encoded by s
"""
def prot(string):
    pattern = {"UUU": "F", "CUU": "L", "AUU": "I", "GUU": "V",
           "UUC": "F", "CUC": "L", "AUC": "I", "GUC": "V",
           "UUA": "L", "CUA": "L", "AUA": "I", "GUA": "V",
           "UUG": "L", "CUG": "L", "AUG": "M", "GUG": "V",
           "UCU": "S", "CCU": "P", "ACU": "T", "GCU": "A",
           "UCC": "S", "CCC": "P", "ACC": "T", "GCC": "A",
           "UCA": "S", "CCA": "P", "ACA": "T", "GCA": "A",
           "UCG": "S", "CCG": "P", "ACG": "T", "GCG": "A",
           "UAU": "Y", "CAU": "H", "AAU": "N", "GAU": "D",
           "UAC": "Y", "CAC": "H", "AAC": "N", "GAC": "D",
           "UAA": "*", "CAA": "Q", "AAA": "K", "GAA": "E",
           "UAG": "*", "CAG": "Q", "AAG": "K", "GAG": "E",
           "UGU": "C", "CGU": "R", "AGU": "S", "GGU": "G",
           "UGC": "C", "CGC": "R", "AGC": "S", "GGC": "G",
           "UGA": "*", "CGA": "R", "AGA": "R", "GGA": "G",
           "UGG": "W", "CGG": "R", "AGG": "R", "GGG": "G"}
    rst = ''.join([pattern[string[i:i+3]] for i in range(0, len(string), 3)])
    return rst[:rst.index("*")]
    
