import random
import matplotlib.pyplot as plt
def random_dna(length):
    return "".join(random.choice("ACGT") for _ in range(length))
genome_len = random.randint(1000, 3000)
genome = random_dna(genome_len)
print(f"Genome length: {genome_len} bases")
num_fragments = 10
min_frag_len = 100
max_frag_len = min(3000, genome_len)

fragments = []         
fragment_lengths = []   

for _ in range(num_fragments):
    L = random.randint(min_frag_len, max_frag_len)
    start = random.randint(0, genome_len - L)
    frag = genome[start:start+L]
    fragments.append(frag)
    fragment_lengths.append(L)

print("Fragment lengths (bp):")
print(fragment_lengths)      

max_len = max(fragment_lengths)
min_len = min(fragment_lengths)


band_positions = [
    0.1 + 0.8 * (max_len - L) / (max_len - min_len + 1e-9)
    for L in fragment_lengths
]

fig, ax = plt.subplots(figsize=(3, 6))

gel_left = 0.4
gel_width = 0.2
gel_bottom = 0.05
gel_height = 0.9
ax.add_patch(
    plt.Rectangle(
        (gel_left, gel_bottom),
        gel_width,
        gel_height,
        color="black"
    )
)

for y in band_positions:
    ax.hlines(y, gel_left + 0.02, gel_left + gel_width - 0.02,
              colors="white", linewidth=3)
for y, L in sorted(zip(band_positions, fragment_lengths), key=lambda x: -x[1]):
    ax.text(gel_left - 0.05, y, f"{L} bp",
            ha="right", va="center", fontsize=8)

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")
plt.tight_layout()
plt.show()