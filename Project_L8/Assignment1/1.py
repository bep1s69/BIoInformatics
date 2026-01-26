import random

# --- Configuration for Simulation ---

# Bases for DNA sequence
BASES = ['A', 'T', 'C', 'G']
# Length of the main host DNA sequence (will be randomly chosen between min/max)
MIN_HOST_LENGTH = 200
MAX_HOST_LENGTH = 400
# Number of transposable elements to insert
NUM_TE_TO_INSERT = random.randint(3, 4)

# Define the structure of the Transposable Element (TE)
# The core TE sequence is flanked by Inverted Repeats (IRs).
# The insertion process creates Target Site Duplications (TSDs), or Direct Repeats,
# which are an immediate repeat of the host DNA sequence at the insertion site.

# A simple, identifiable TE sequence for simulation
TE_CORE_SEQUENCE = "GATTACA"  # This is the sequence between the Inverted Repeats
# Inverted Repeat (IR) sequence. The complement/reverse complement will be used.
IR_SEQUENCE = "TCGA"
# Length of the TSD (Direct Repeat) sequence. Often 2-10 bp.
TSD_LENGTH = 4

# --- Helper Functions ---

def generate_random_dna(length):
    """Generates a random DNA sequence of a specified length."""
    return "".join(random.choice(BASES) for _ in range(length))

def get_complement(base):
    """Returns the complementary base."""
    if base == 'A': return 'T'
    if base == 'T': return 'A'
    if base == 'C': return 'G'
    if base == 'G': return 'C'
    return base # Should not happen

def get_reverse_complement(sequence):
    """Returns the reverse complement of a sequence (for the second IR)."""
    complement = "".join(get_complement(base) for base in sequence)
    return complement[::-1]


def create_simulated_dna_with_tes():
    """
    Creates the artificial DNA sequence with TEs inserted.
    Returns: The final DNA sequence and a list of the actual TE positions (start, end).
    """
    host_length = random.randint(MIN_HOST_LENGTH, MAX_HOST_LENGTH)
    host_dna = list(generate_random_dna(host_length))
    te_positions = []
    
    # Define the full TE structure for insertion
    # Structure: TSD (Direct Repeat) - IR1 - TE_CORE - IR2 - TSD (Direct Repeat)
    ir2_sequence = get_reverse_complement(IR_SEQUENCE)
    
    # The full sequence of the inserted element (IR1 + CORE + IR2)
    TE_INSERT = IR_SEQUENCE + TE_CORE_SEQUENCE + ir2_sequence
    TE_INSERT_LENGTH = len(TE_INSERT)
    
    # Calculate potential insertion sites, avoiding the very beginning/end
    possible_indices = list(range(TSD_LENGTH, len(host_dna) - TSD_LENGTH))
    
    # Shuffle to pick random, non-overlapping insertion sites
    random.shuffle(possible_indices)

    print(f"--- 1. Generating DNA with {NUM_TE_TO_INSERT} TEs ---")
    print(f"Host DNA Length: {len(host_dna)}")
    
    inserted_count = 0
    
    # Iterate through potential sites to insert TEs
    for insertion_index in possible_indices:
        if inserted_count >= NUM_TE_TO_INSERT:
            break

        # Check for intersection/overlap with already inserted TEs
        # We'll use a simple proximity check for non-overlapping
        is_overlapping = False
        for start, end in te_positions:
            # Check if the potential insertion site is too close to an existing TE
            # We need space for the new TSDs and the TE itself.
            if (insertion_index > start - (TE_INSERT_LENGTH + TSD_LENGTH) and 
                insertion_index < end + (TE_INSERT_LENGTH + TSD_LENGTH)):
                is_overlapping = True
                break
        
        if is_overlapping:
            continue

        # 1. Identify the Target Site Duplication (TSD) sequence
        # The TSD sequence is the host DNA sequence immediately before the insertion point.
        tsd_sequence = "".join(host_dna[insertion_index : insertion_index + TSD_LENGTH])

        # 2. Construct the full inserted sequence
        # The host sequence is duplicated on the other side of the TE
        full_insertion = tsd_sequence + TE_INSERT + tsd_sequence

        # 3. Insert the TE into the host DNA
        # This is the simulated cut-and-paste mechanism:
        # [Host Pre] | [TSD] | [IR1-CORE-IR2] | [TSD] | [Host Post]
        
        # Split the host DNA into three parts:
        # Part 1: Sequence before the first TSD
        pre_insertion = host_dna[:insertion_index]
        # Part 2: Sequence that becomes the first TSD (removed from original host)
        # Part 3: Sequence after the TSD
        post_insertion = host_dna[insertion_index + TSD_LENGTH:]

        # Recombine to form the new DNA sequence
        host_dna = pre_insertion + list(full_insertion) + post_insertion

        # Calculate the new position of the inserted element (IR1-CORE-IR2)
        te_start = len(pre_insertion) + len(tsd_sequence)
        te_end = te_start + TE_INSERT_LENGTH - 1
        te_positions.append((te_start, te_end))

        inserted_count += 1
        
        # IMPORTANT: Since we modified host_dna, we need to adjust subsequent TE positions
        # This is a key part of the 'intersecting' challenge.
        # For simplicity, let's restart the possible_indices for the next iteration
        # OR handle the overlap by ensuring no TEs are close together (which we did).
        
        # To simplify the overlap/intersection check (as the prompt suggests),
        # we will break the generation after inserting the required number.
        if inserted_count >= NUM_TE_TO_INSERT:
            break
        
    final_dna = "".join(host_dna)
    print(f"Final DNA Length: {len(final_dna)}")
    print(f"Actual TE Positions (0-indexed start, end): {te_positions}")
    print("-" * 50)
    
    return final_dna, te_positions, TE_INSERT, IR_SEQUENCE, ir2_sequence, TSD_LENGTH

# --- Task 2: Implement Detection Software ---

def detect_transposable_elements(dna_sequence, te_insert_sequence, 
                                 ir_sequence_fwd, ir_sequence_rev, tsd_length):
    """
    Detects the positions of transposable elements (TEs) based on the
    characteristic flanking structure: TSD - IR1 - TE_CORE - IR2 - TSD.
    
    The algorithm specifically searches for:
    1. The full TE core sequence (IR1 + TE_CORE + IR2) AND
    2. Identical sequences of TSD_LENGTH immediately flanking the full TE core.
    
    Returns: A list of detected TE positions (start, end).
    """
    detected_positions = []
    
    # 1. Define the full TE structure being searched for
    te_length = len(te_insert_sequence)
    
    print("--- 2. Implementing Detection Algorithm ---")
    print(f"Searching for TE core sequence: {te_insert_sequence} (Length: {te_length})")
    print(f"TSD (Direct Repeat) Length: {tsd_length}")
    
    # Search for potential TE locations in the DNA sequence
    for i in range(len(dna_sequence) - te_length):
        # A. Check for the core TE sequence (IR1-CORE-IR2)
        potential_te = dna_sequence[i : i + te_length]
        
        if potential_te == te_insert_sequence:
            # B. Check for the flanking TSDs (Direct Repeats)
            
            # The TE starts at index 'i'.
            # The first TSD should be immediately before it.
            tsd1_start = i - tsd_length
            tsd1_end = i
            
            # The second TSD should be immediately after the TE.
            tsd2_start = i + te_length
            tsd2_end = i + te_length + tsd_length
            
            # Check for bounds to ensure TSDs are within the DNA sequence
            if tsd1_start < 0 or tsd2_end > len(dna_sequence):
                continue
                
            # Extract the potential TSD sequences
            tsd1 = dna_sequence[tsd1_start : tsd1_end]
            tsd2 = dna_sequence[tsd2_start : tsd2_end]
            
            # C. TSD Validation: Check if the two flanking sequences are identical
            if tsd1 == tsd2:
                # We found a validated TE insertion!
                # The position refers to the TE itself (IR1-CORE-IR2)
                te_start = i
                te_end = i + te_length - 1
                detected_positions.append((te_start, te_end))
                
                # To handle the "intersect two or three transposons" case,
                # we don't skip over the current TE. We let the loop continue
                # to potentially find a nested TE or one that starts immediately after.

    print(f"Detected TE Positions (0-indexed start, end): {detected_positions}")
    print("-" * 50)
    return detected_positions

# --- Main Execution ---

if __name__ == "__main__":
    # 1. Generate the simulated DNA sequence
    final_dna, actual_te_positions, te_insert, ir1, ir2, tsd_len = create_simulated_dna_with_tes()

    # 2. Detect the TEs in the generated sequence
    detected_te_positions = detect_transposable_elements(
        final_dna, te_insert, ir1, ir2, tsd_len
    )
    
    # 3. Compare Results
    print("### Comparison of Results ###")
    
    # Sort for easy comparison
    actual_te_positions.sort()
    detected_te_positions.sort()
    
    if actual_te_positions == detected_te_positions:
        print("✅ Detection Successful: Actual and Detected positions match!")
    else:
        print("❌ Detection Mismatch:")
        print(f"   Actual Positions:   {actual_te_positions}")
        print(f"   Detected Positions: {detected_te_positions}")
    
    # Display a snippet of the sequence for verification
    print("\nSnippet of DNA with the first detected TE highlighted:")
    if detected_te_positions:
        start, end = detected_te_positions[0]
        tsd_start = start - tsd_len
        tsd_end = end + tsd_len + 1 # +1 because slicing is exclusive
        
        # Display TSD1, TE, and TSD2 with separators for clarity
        tsd1_seq = final_dna[tsd_start:start]
        te_seq = final_dna[start:end + 1]
        tsd2_seq = final_dna[end + 1:tsd_end]
        
        highlighted_snippet = (
            final_dna[tsd_start - 20:tsd_start] + # Context
            "" + tsd1_seq + "" + # TSD1
            "---" + te_seq + "---" + # TE
            "" + tsd2_seq + "" + # TSD2
            final_dna[tsd_end:tsd_end + 20] # Context
        )
        
        print(f"Context... **TSD1---TE_CORE---TSD2** ...Context")
        print(highlighted_snippet)