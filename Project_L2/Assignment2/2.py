S = "TACGTGGGGGGGAGCTATCTACTGACTTACGACTAGCTAGCTGCATCATCGATCGA"

# Function to extract dinucleotides
def get_dinucleotides(sequence):
    # Use a sliding window of size 2
    dinucleotides = {sequence[i:i+2] for i in range(len(sequence) - 1)}
    return sorted(dinucleotides)

# Function to extract trinucleotides
def get_trinucleotides(sequence):
    # Use a sliding window of size 3
    trinucleotides = {sequence[i:i+3] for i in range(len(sequence) - 2)}
    return sorted(trinucleotides)

# Get and print the results
dinucleotides = get_dinucleotides(S)
trinucleotides = get_trinucleotides(S)

print("Dinucleotides:", ", ".join(dinucleotides))
print("Trinucleotides:", ", ".join(trinucleotides))
