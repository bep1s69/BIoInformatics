from Bio import Entrez, SeqIO
import matplotlib.pyplot as plt
 
accessions = [
    "NC_007366",  # example: Influenza A virus (H3N2) HA segment
    "NC_007362",  # Influenza A H5N1 HA segment
    "NC_002018",  # Influenza A PR8 segment
    "NC_007367",
    "NC_007368",
    "NC_007369",
    "NC_007370",
    "NC_007371",
    "NC_007372",
    "NC_007373",
]

# Enzymes from assignment 1: (recognition_sequence, cut_pos_inside_motif)
enzymes = {
    "EcoRI":  ("GAATTC", 1),   # G^AATTC
    "BamHI":  ("GGATCC", 1),   # G^GATCC
    "HindIII":("AAGCTT", 1),   # A^AGCTT
    "TaqI":   ("TCGA",   1),   # T^CGA
    "HaeIII": ("GGCC",   2)    # GG^CC
}

TOL = 5  # bp tolerance when deciding “same band” across genomes

# HELPERS: download & digest
def fetch_sequence(accession):
    """Download one nucleotide sequence from NCBI."""
    handle = Entrez.efetch(db="nucleotide", id=accession,
                           rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    seq = str(record.seq).upper()
    return record.description, seq


def find_cuts(seq, recog, cutpos):
    """Return list of cut positions (0-based indices in full sequence)."""
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
    """From cut positions produce fragment length list."""
    cuts = sorted([0] + cuts + [seq_length])
    return [cuts[i+1] - cuts[i] for i in range(len(cuts) - 1)]


def digest_with_all_enzymes(seq):
    """Union of all cut sites for the 5 enzymes → one pattern per genome."""
    all_cuts = []
    for recog, cutpos in enzymes.values():
        all_cuts.extend(find_cuts(seq, recog, cutpos))
    all_cuts = sorted(set(all_cuts))
    return compute_fragments(len(seq), all_cuts)

#  DOWNLOAD 10 GENOMES & DIGEST

genome_fragments = {}   # accession → fragment list
genome_descriptions = {}

for acc in accessions:
    print(f"Downloading {acc} ...")
    desc, seq = fetch_sequence(acc)
    genome_descriptions[acc] = desc
    frags = digest_with_all_enzymes(seq)
    genome_fragments[acc] = frags

    print("  Description:", desc)
    print("  Length:", len(seq), "bp")
    print("  Fragments after all 5 enzymes:", frags, "\n")

#  GEL FOR EACH GENOME
def draw_gel(digests, title):
    """
    digests: dict {label: [fragment_sizes]}
    One lane per label.
    """
    if not digests:
        print("No data to plot.")
        return

    max_frag = max(max(v) for v in digests.values())
    y_max = max_frag + 200

    plt.figure(figsize=(10, 8))
    lane = 1
    for label, frags in digests.items():
        for f in sorted(frags, reverse=True):
            y = y_max - f
            plt.hlines(y, lane - 0.3, lane + 0.3, linewidth=6)
        plt.text(lane, y_max + 80, label, ha="center",
                 va="bottom", rotation=90, fontsize=8)
        lane += 1

    plt.gca().invert_yaxis()
    plt.xticks([])
    plt.ylabel("Migration distance (shorter fragments → lower)")
    plt.title(title)
    plt.tight_layout()
    plt.show()


# Original gel: 10 genomes, all bands
draw_gel(genome_fragments, "a) All genomes – all fragments (5-enzyme digest)")

#  REMOVE BANDS COMMON TO ALL
def is_common_band(size, all_frag_lists, tol=TOL):
    """
    True if this fragment size appears in ALL fragment lists
    within ±tol bp.
    """
    for frags in all_frag_lists:
        if not any(abs(size - f) <= tol for f in frags):
            return False
    return True


all_frag_lists = list(genome_fragments.values())
difference_fragments = {}

for acc, frags in genome_fragments.items():
    unique_frags = [
        f for f in frags
        if not is_common_band(f, all_frag_lists, tol=TOL)
    ]
    difference_fragments[acc] = unique_frags
    print(f"{acc} unique fragments (after removing common bands):")
    print(unique_frags, "\n")

draw_gel(difference_fragments,
         "c) General gel – only bands that differ between genomes")