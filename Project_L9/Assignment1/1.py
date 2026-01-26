from Bio import Entrez, SeqIO
import matplotlib.pyplot as plt
 

accession = "NM_001101"

handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
record = SeqIO.read(handle, "fasta")
handle.close()

sequence = str(record.seq).upper()
seq_len = len(sequence)

print(f"Downloaded: {record.id}")
print(f"Description: {record.description}")
print(f"Sequence length: {seq_len} bp\n")

# Restriction enzymes
enzymes = {
    "EcoRI":  ("GAATTC", 1),
    "BamHI":  ("GGATCC", 1),
    "HindIII":("AAGCTT", 1),
    "TaqI":   ("TCGA",   1),
    "HaeIII": ("GGCC",   2)
}

def find_cuts(seq, recog, cutpos):
    cuts = []
    start = 0
    while True:
        idx = seq.find(recog, start)
        if idx == -1:
            break
        cuts.append(idx + cutpos)
        start = idx + 1
    return cuts

def compute_fragments(seq_length, cuts):
    cuts = sorted([0] + cuts + [seq_length])
    return [cuts[i+1] - cuts[i] for i in range(len(cuts)-1)]

# Compute digests
all_digests = {}

for name, (recog, cutpos) in enzymes.items():
    cuts = find_cuts(sequence, recog, cutpos)
    fragments = compute_fragments(seq_len, cuts)
    all_digests[name] = fragments

    print(f"=== {name} ===")
    print("Recognition:", recog)
    print("Cut positions:", [c+1 for c in cuts])  # 1-based
    print("Fragment lengths:", fragments, "\n")

# Multi-digest (all enzymes together)
all_cuts = []
for name, (recog, cutpos) in enzymes.items():
    all_cuts += find_cuts(sequence, recog, cutpos)
all_cuts = sorted(set(all_cuts))
multi = compute_fragments(seq_len, all_cuts)
all_digests["All_5"] = multi

print("=== All 5 Enzymes Digest ===")
print("Cut positions:", [c+1 for c in all_cuts])
print("Fragment lengths:", multi)

#  Gel simulation
def draw_gel(digests):
    max_frag = max(max(f) for f in digests.values())
    y_max = max_frag + 200

    plt.figure(figsize=(7, 8))
    lane = 1
    for enzyme, fragments in digests.items():
        for f in sorted(fragments, reverse=True):
            y = y_max - f
            plt.hlines(y, lane - 0.3, lane + 0.3, linewidth=8)
        plt.text(lane, y_max + 80, enzyme, ha='center', rotation=45)
        lane += 1

    plt.gca().invert_yaxis()
    plt.xticks([])
    plt.ylabel("Migration (shorter fragments → further)")
    plt.title("Simulated Restriction Digest Gel")
    plt.tight_layout()
    plt.show()


draw_gel(all_digests)