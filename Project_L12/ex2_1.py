import json
import random
import os

def load_json(filename):
    """Loads JSON data from a file."""
    if not os.path.exists(filename):
        print(f"Error: {filename} not found. Please run previous assignments first.")
        return None
    with open(filename, 'r') as f:
        return json.load(f)

def synthesize_sequence(matrix, start_state, length=20):
    """
    Generates a sequence by walking through the Markov Chain.
    """
    current_state = start_state
    sequence = [current_state]

    for _ in range(length - 1):
        if current_state not in matrix:
            break # Stop if we hit a dead end
            
        transitions = matrix[current_state]
        possible_next_states = list(transitions.keys())
        probabilities = list(transitions.values())

        if not possible_next_states:
            break

        # Weighted random selection based on probabilities
        next_state = random.choices(possible_next_states, weights=probabilities, k=1)[0]
        sequence.append(next_state)
        current_state = next_state
        
    return sequence

def generate_dna_outputs(count=3):
    """Generates multiple DNA sequences."""
    data = load_json("dna_transitions.json")
    if not data: return []

    results = []
    keys = list(data.keys())
    
    for i in range(count):
        start_node = random.choice(keys)
        # Generate 50 nucleotides
        seq_list = synthesize_sequence(data, start_node, length=50)
        dna_str = "".join(seq_list)
        results.append(dna_str)
        print(f"[DNA {i+1}] {dna_str[:30]}...") # Print preview
        
    return results

def generate_text_outputs(count=3):
    """Generates multiple Text sequences."""
    data = load_json("word_transitions.json")
    if not data: return []

    legend = data["legend"]
    matrix = data["matrix"]
    results = []
    keys = list(matrix.keys())
    
    for i in range(count):
        start_symbol = random.choice(keys)
        # Generate 20 words
        symbol_seq = synthesize_sequence(matrix, start_symbol, length=20)
        
        # Decode symbols back to words
        word_seq = [legend[sym] for sym in symbol_seq if sym in legend]
        text_str = " ".join(word_seq)
        results.append(text_str)
        print(f"[Text {i+1}] {text_str[:40]}...") # Print preview

    return results

if __name__ == "__main__":
    print("--- Synthesizing Sequences ---")
    
    # Generate 3 outputs for each type
    dna_results = generate_dna_outputs(3)
    text_results = generate_text_outputs(3)
    
    # Structure the output data
    output_data = {
        "DNA_Sequences": dna_results,
        "Text_Sequences": text_results
    }
    
    # Save to JSON
    filename = "synthesized_sequences.json"
    with open(filename, 'w') as f:
        json.dump(output_data, f, indent=4)
        
    print(f"\nSaved all outputs to: {os.path.abspath(filename)}")