from Bio import pairwise2
from Bio.Align import substitution_matrices

matrix=substitution_matrices.load('BLOSUM45')

gapopen=-1
gapextend=-1

# seqA='AGCTTGCAGTTACGGCTAGCT'
# seqB='AGCTTGCAGTTACAGCTAGCT'
def pairwiseseqalign(seqA,seqB):
   aln=pairwise2.align.globalds(seqA,seqB,matrix,gapopen,gapextend)
   stringer=str((pairwise2.format_alignment(*aln[0])))
   identify=str(round(((stringer.count ('|'))/(aln[0].end))*100,2))
   print(identify)
   print(stringer)
   print((pairwise2.format_alignment(*aln[0])))
   print(aln[0].end)
#    print(aln[1].end)

pairwiseseqalign('AGCTTGCAGTTACGGCTAGCT','AGCTTGCAGTTACAGCTAGCT')
pairwiseseqalign('AGTACCGTACCGTACCGTACCGTACCGTACCGTACCGTACCGTACCGTA','TACGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA')


