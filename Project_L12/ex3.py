import string
import json
import os

def get_sample_text():
    """
    Returns a random English text of approximately 300 characters.
    """
    text = (
        "Bioinformatics is an interdisciplinary field that develops methods and "
        "software tools for understanding biological data, in particular when the "
        "data sets are large and complex. As an interdisciplinary field of science, "
        "bioinformatics combines biology, computer science, information engineering, "
        "mathematics and statistics to analyze and interpret the biological data."
    )
    return text

def preprocess_and_tokenize(text):
    """
    Converts text to lowercase, removes punctuation, and splits into words.
    """
    # Remove punctuation using a translation table
    translator = str.maketrans('', '', string.punctuation)
    clean_text = text.translate(translator).lower()
    
    # Split by whitespace into a list of words
    words = clean_text.split()
    return words

def map_words_to_symbols(words):
    """
    Assigns a unique ASCII symbol to each unique word.
    Returns:
        word_to_sym (dict): maps 'biology' -> 'A'
        sym_to_word (dict): maps 'A' -> 'biology' (for the legend)
    """
    unique_words = sorted(list(set(words)))
    
    # We use printable ASCII characters starting from '!' (ASCII 33)
    # to avoid control characters.
    word_to_sym = {}
    sym_to_word = {}
    
    for i, word in enumerate(unique_words):
        # chr(33 + i) gives us symbols like !, ", #, $, %, etc.
        symbol = chr(33 + i) 
        word_to_sym[word] = symbol
        sym_to_word[symbol] = word
        
    return word_to_sym, sym_to_word

def compute_word_transition_matrix(words, word_to_sym):
    """
    Computes probabilities of Symbol(Word A) -> Symbol(Word B).
    """
    # Initialize counts
    unique_syms = list(word_to_sym.values())
    counts = {s1: {s2: 0 for s2 in unique_syms} for s1 in unique_syms}

    # Count transitions
    for i in range(len(words) - 1):
        curr_sym = word_to_sym[words[i]]
        next_sym = word_to_sym[words[i+1]]
        counts[curr_sym][next_sym] += 1

    # Normalize to probabilities
    matrix = {}
    for s1 in unique_syms:
        total = sum(counts[s1].values())
        matrix[s1] = {}
        for s2 in unique_syms:
            if total > 0 and counts[s1][s2] > 0:
                # Only storing non-zero probabilities to keep JSON clean
                matrix[s1][s2] = round(counts[s1][s2] / total, 4)
    
    return matrix

# --- Main Execution ---
if __name__ == "__main__":
    # 1. Get Text
    text = get_sample_text()
    print(f"Input Text ({len(text)} chars):\n{text}\n")

    # 2. Tokenize
    word_list = preprocess_and_tokenize(text)
    print(f"Total words found: {len(word_list)}")

    # 3. Create Symbol Mapping
    w2s, s2w = map_words_to_symbols(word_list)
    
    # 4. Compute Matrix
    transition_matrix = compute_word_transition_matrix(word_list, w2s)

    # 5. Save to JSON (Including the Legend!)
    output_data = {
        "legend": s2w,          # Helps you decode which symbol is which word
        "matrix": transition_matrix
    }

    filename = "word_transitions.json"
    with open(filename, 'w') as f:
        json.dump(output_data, f, indent=4)
        
    print(f"File saved: {os.path.abspath(filename)}")
    print("\nSample of the Legend (Symbol -> Word):")
    sample_items = list(s2w.items())[:5]
    for sym, word in sample_items:
        print(f"  '{sym}' : {word}")