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

if __name__ == "__main__":
    with open("../data/rosalind_subs.txt", 'r') as f:
        string1 = f.readline().strip()
        string2 = f.readline().strip()
    loc = subs(string1, string2)
    for i in loc:
        print(i, end=" ")



"""
Calculating Protein Mass
URL: http://rosalind.info/problems/prtm/

Given: A protein string P of length at most 1000 aa.
Return: The total weight of P. Consult the monoisotopic mass table.
"""



"""
RNA Splicing
URL: http://rosalind.info/problems/splc/

Given: A DNA string s (of length at most 1 kbp) and a collection of substrings of s acting as introns. All strings are given in FASTA format.
Return: A protein string resulting from transcribing and translating the exons of s. (Note: Only one solution will exist for the dataset provided.)
"""



"""
Locating Restriction Sites
URL: http://rosalind.info/problems/revp/

Given: A DNA string of length at most 1 kbp in FASTA format.
Return: The position and length of every reverse palindrome in the string having length between 4 and 12. You may return these pairs in any order.
"""



"""
Transitions and Transversions
URL: http://rosalind.info/problems/tran/

Given: Two DNA strings s1 and s2 of equal length (at most 1 kbp).
Return: The transition/transversion ratio R(s1,s2).
"""



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
