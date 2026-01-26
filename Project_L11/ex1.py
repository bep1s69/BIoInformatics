import math
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

S1_island = "ATCGATTCGATATCATACACGTAT"          # (+) Model
S2_non_island = "CTCGACTAGTATGAAGTCCACGCTTG"    # (-) Model
S_new = "CAGGTTGGAAACGTAA"                      # Sequence to test

bases = ['A', 'C', 'G', 'T']

#HELPER FUNCTIONS HERE 

def get_transition_counts(sequence):
    """
    Initializes a 4x4 matrix and counts transitions (e.g., A->C).
    """
    counts = {b1: {b2: 0 for b2 in bases} for b1 in bases}
    for i in range(len(sequence) - 1):
        c1 = sequence[i]
        c2 = sequence[i+1]
        counts[c1][c2] += 1
    return counts

def calculate_probabilities(counts):
    """
    Converts counts to probabilities using Laplace smoothing (+1).
    """
    probs = {b1: {b2: 0.0 for b2 in bases} for b1 in bases}
    for b1 in bases:
        # Sum of row + 4 (for the +1 to each of the 4 bases)
        row_sum = sum(counts[b1].values()) + 4 
        for b2 in bases:
            probs[b1][b2] = (counts[b1][b2] + 1) / row_sum
    return probs

#CALCULATE THE MATRIX

# Count transitions
counts_plus = get_transition_counts(S1_island)
counts_minus = get_transition_counts(S2_non_island)

# Convert to probabilities
probs_plus = calculate_probabilities(counts_plus)
probs_minus = calculate_probabilities(counts_minus)

# Calculate Log-Likelihoods
ll_data = []
for b1 in bases:
    row = []
    for b2 in bases:
        p_plus = probs_plus[b1][b2]
        p_minus = probs_minus[b1][b2]
        score = math.log2(p_plus / p_minus)
        row.append(score)
    ll_data.append(row)

# Create a Pandas DataFrame for easier handling and plotting
ll_df = pd.DataFrame(ll_data, index=bases, columns=bases)

print("--- Calculated Log-Likelihood Matrix ---")
print(ll_df.round(3))
print("\n")

#PLOTS

print("Generating plot...")

plt.figure(figsize=(8, 6))

# Create the heatmap
sns.heatmap(ll_df, annot=True, fmt=".3f", cmap="coolwarm", center=0,
            linewidths=1, linecolor='black', cbar_kws={'label': 'Log-Likelihood Score'})

plt.title("Log-Likelihood Matrix\n(CpG Island vs Background)", fontsize=14)
plt.xlabel("To Nucleotide", fontsize=12)
plt.ylabel("From Nucleotide", fontsize=12)

# Save the plot to the current folder
output_filename = "log_likelihood_heatmap.png"
plt.savefig(output_filename, dpi=300, bbox_inches='tight')

# Show the plot (optional, depends on your environment)
# plt.show() 

print(f"Plot saved successfully as '{output_filename}'\n")
total_score = 0
print(f"Scoring Sequence S: {S_new}")
print("-" * 30)

for i in range(len(S_new) - 1):
    curr_n = S_new[i]
    next_n = S_new[i+1]
    
    # Look up score in our DataFrame
    step_score = ll_df.loc[curr_n, next_n]
    total_score += step_score
    
    print(f"{curr_n} -> {next_n} : {step_score: .4f}")

print("-" * 30)
print(f"Total Log-Likelihood Score: {total_score:.4f}\n")

#conclusion

if total_score > 0:
    print(">>> Result: POSITIVE (+)")
    print(">>> The sequence belongs to a CpG Island.")
else:
    print(">>> Result: NEGATIVE (-)")
    print(">>> The sequence does NOT belong to a CpG Island.")