from collections import defaultdict

def find_repeats(dna_sequence, min_len=3, max_len=6, min_repeats=2):
    """
    Finds repeating subsequences between min_len and max_len (inclusive)
    that repeat at least min_repeats times.
    Returns a dictionary: subseq -> list of positions (0-based).
    """
    seq_len = len(dna_sequence)
    repeat_dict = defaultdict(list)

    for length in range(min_len, max_len+1):
        seen = {}
        for start in range(seq_len - length + 1):
            subseq = dna_sequence[start:start+length]
            repeat_dict[subseq].append(start)

    # Filter only those with at least min_repeats occurrences
    filtered = {seq: pos for seq, pos in repeat_dict.items() if len(pos) >= min_repeats}
    return filtered

def main():
    
    dna_sequence = (
        "TGCATGCTAGCTGTCATGCTTAGCTAGTGCATGCTAGCTGTCATGCTTAAAGCTAGCT"
        "GTCATGCTT" * 5
    )
    print('Sequence length:', len(dna_sequence))

    repeats = find_repeats(dna_sequence, min_len=3, max_len=6, min_repeats=2)
    print("\nDetected repeating sequences (3b-6b, at least 2 repeats):\n")

    for subseq, positions in sorted(repeats.items(),
                                    key=lambda x: (-len(x[1]), x[0])):
        print(f"Subsequence: {subseq}, Count: {len(positions)}, Positions: {positions}")

if __name__ == "__main__":
    main()
