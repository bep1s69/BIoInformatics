import itertools
from collections import Counter

S = "TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"
nuc = 'ACGT'

# Dinucleotides
all_di = [''.join(p) for p in itertools.product(nuc, repeat=2)]
di_counts = Counter(S[i:i+2] for i in range(len(S)-1))
total_di = len(S) - 10
print("Dinucleotides:")
for combo in sorted(all_di):
    cnt = di_counts.get(combo, 0)
    perc = (cnt / total_di) * 100 if total_di > 0 else 0
    print(f"{combo}: {perc:.2f}%")

# Trinucleotides
all_tri = [''.join(p) for p in itertools.product(nuc, repeat=3)]
tri_counts = Counter(S[i:i+3] for i in range(len(S)-2))
total_tri = len(S) - 2
print("\nTrinucleotides:")
for combo in sorted(all_tri):
    cnt = tri_counts.get(combo, 0)
    perc = (cnt / total_tri) * 100 if total_tri > 0 else 0
    print(f"{combo}: {perc:.2f}%")