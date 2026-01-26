import random
import json
import numpy as np
import os

def generate_dna_sequence(length=50):
    """Generates a random DNA sequence of given length."""
    return "".join(random.choice("ACGT") for _ in range(length))

def compute_transition_matrix(sequence):
    """
    Computes the transition probability matrix for a DNA sequence.
    Returns a dictionary structure suitable for JSON export.
    """
    nucleotides = ['A', 'C', 'G', 'T']
    
    # Initialize a 4x4 matrix with zeros
    # Rows = Current State, Columns = Next State
    counts = {n1: {n2: 0 for n2 in nucleotides} for n1 in nucleotides}
    
    # Count transitions (loop up to the second to last character)
    for i in range(len(sequence) - 1):
        current_n = sequence[i]
        next_n = sequence[i+1]
        counts[current_n][next_n] += 1

    # Convert counts to probabilities (Normalize)
    transition_matrix = {}
    
    for start_n in nucleotides:
        total_transitions = sum(counts[start_n].values())
        transition_matrix[start_n] = {}
        
        for end_n in nucleotides:
            if total_transitions > 0:
                prob = counts[start_n][end_n] / total_transitions
            else:
                prob = 0.0 # Handle case where a nucleotide never appears
            
            # Rounding for cleaner JSON output
            transition_matrix[start_n][end_n] = round(prob, 4)

    return transition_matrix

def save_to_json(data, filename="dna_transitions.json"):
    """Saves the dictionary to a JSON file."""
    with open(filename, 'w') as f:
        json.dump(data, f, indent=4)
    print(f"File saved successfully: {os.path.abspath(filename)}")

# --- Main Execution ---
if __name__ == "__main__":
    # 1. Generate Sequence
    dna_seq = generate_dna_sequence(50)
    print(f"Generated Sequence:\n{dna_seq}\n")

    # 2. Compute Matrix
    matrix_data = compute_transition_matrix(dna_seq)

    # 3. Print Matrix to Console
    print("Transition Matrix (Probability of Row -> Col):")
    print(f"{'':<4} {'A':<6} {'C':<6} {'G':<6} {'T':<6}")
    for start_n, transitions in matrix_data.items():
        print(f"{start_n:<4}", end="")
        for end_n in ['A', 'C', 'G', 'T']:
            print(f"{transitions[end_n]:<6.4f} ", end="")
        print() # Newline

    # 4. Save to JSON
    save_to_json(matrix_data)