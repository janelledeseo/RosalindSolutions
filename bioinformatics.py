### ICA Session 01: 1/23/25
## FIRST ATTEMPT
## For a given 'dna.sequence', write one basic function and calculate GC content of a string

# To count the occurrences of G and C
g_count = dna_sequence.count('G')
c_count = dna_sequence.count('C')
    
# To calculate the GC content
gc_content = (g_count + c_count) / len(dna_sequence) * 100
gc_content

### After in-class review:
## EXAMPLE:
sequence1 = 'AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC'
def gc_content(seq):
  return (seq.count('G') + seq.count('C'))/len(seq)
