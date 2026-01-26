import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re

# ---------------------------------------------------------
# 1. DATA PREPARATION
# ---------------------------------------------------------

text_eminescu = """
A fost odata ca-n povesti a fost ca niciodata din rude mari imparatesti o prea frumoasa fata.
Si era una la parinti si mandra-n toate cele cum e fecioara intre sfinti si luna intre stele.
Cobori in jos luceafar bland alunecand pe o raza patrunde-n casa si in gand si viata-mi lumineaza.
Mai am un singur dor in linistea de seara sa ma lasati sa mor la marginea de mare.
"""

text_stanescu = """
Leoaica tanara iubirea mi-ai sarit in fata. Ma pandise-n incordare mai demult.
Coltii albi mi i-a infipt in fata m-a muscat leoaica azi de fata.
Si de-aude ploaia plangand la geamuri si vantul cum suiera prin ramuri.
E un fel de a fi auriu si un fel de a fi curat.
Pasul tau de domnisoara, pasul tau de pasare, pasul tau de toamna pe trotuare.
"""

# The "ACCUSED" Text (Mix: Eminescu -> Stanescu -> Eminescu)
text_accused = (
    "cobori in jos luceafar bland alunecand pe o raza " \
    "leoaica tanara iubirea mi-ai sarit in fata " \
    "patrunde-n casa si in gand si viata-mi lumineaza"
)

valid_chars = "abcdefghijklmnopqrstuvwxyz "

# ---------------------------------------------------------
# 2. PRE-PROCESSING & TRAINING
# ---------------------------------------------------------

def clean_text(text):
    text = text.lower()
    text = re.sub(r'[^a-z ]', '', text)
    text = re.sub(r'\s+', ' ', text).strip()
    return text

clean_eminescu = clean_text(text_eminescu)
clean_stanescu = clean_text(text_stanescu)
clean_accused = clean_text(text_accused)

def get_markov_matrix(text, chars):
    counts = {c1: {c2: 1 for c2 in chars} for c1 in chars} # Laplace +1
    for i in range(len(text) - 1):
        counts[text[i]][text[i+1]] += 1
            
    probs = {}
    for c1 in chars:
        total = sum(counts[c1].values())
        probs[c1] = {c2: counts[c1][c2] / total for c2 in chars}
    return probs

# Train models
prob_eminescu = get_markov_matrix(clean_eminescu, valid_chars)
prob_stanescu = get_markov_matrix(clean_stanescu, valid_chars)

# Build Log-Likelihood Matrix
ll_matrix = {c1: {} for c1 in valid_chars}
for c1 in valid_chars:
    for c2 in valid_chars:
        p_e = prob_eminescu[c1][c2]
        p_s = prob_stanescu[c1][c2]
        ll_matrix[c1][c2] = math.log2(p_e / p_s)

ll_df = pd.DataFrame(ll_matrix).T

# ---------------------------------------------------------
# 3. SLIDING WINDOW ANALYSIS (Smoother)
# ---------------------------------------------------------

# INCREASED WINDOW SIZE FOR SMOOTHING
window_size = 30 
scores = []
positions = []

print(f"Scanning accused text (Length: {len(clean_accused)} chars)...")
print(f"Window Size: {window_size} (Increased for smoothing)\n")

for i in range(len(clean_accused) - window_size + 1):
    window_seq = clean_accused[i : i + window_size]
    window_score = 0
    for j in range(len(window_seq) - 1):
        c1 = window_seq[j]
        c2 = window_seq[j+1]
        window_score += ll_df.loc[c1, c2]
    scores.append(window_score)
    # Record the middle position of the window for plotting
    positions.append(i + window_size // 2)

# Convert to numpy arrays for easier plotting logic
scores_np = np.array(scores)
positions_np = np.array(positions)

# ---------------------------------------------------------
# 4. VISUALIZATION (Smooth Line Plot with Fill)
# ---------------------------------------------------------

plt.figure(figsize=(12, 6))

# Plot the main score line
plt.plot(positions_np, scores_np, color='black', linewidth=1.5, alpha=0.6)

# Fill area above 0 (Eminescu) with Red
plt.fill_between(positions_np, scores_np, 0, 
                 where=(scores_np >= 0), interpolate=True, color='red', alpha=0.5, label='Eminescu Style (+)')

# Fill area below 0 (Stanescu) with Blue
plt.fill_between(positions_np, scores_np, 0, 
                 where=(scores_np <= 0), interpolate=True, color='blue', alpha=0.5, label='Stănescu Style (-)')

# Add baseline and labels
plt.axhline(0, color='black', linewidth=1, linestyle='--')
plt.title(f"Plagiarism Detection: Smoothed Sliding Window ({window_size} chars)", fontsize=14)
plt.xlabel("Approx. Character Position in Text", fontsize=12)
plt.ylabel("Log-Likelihood Score (Smoothed)", fontsize=12)

# Clean legend instead of overlapping text
plt.legend(loc='best', frameon=True)

plt.tight_layout()
plt.savefig("plagiarism_detection_smooth.png")
print("Saved smooth analysis plot to 'plagiarism_detection_smooth.png'")
plt.show()