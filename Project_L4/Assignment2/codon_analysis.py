import re
from collections import Counter
import matplotlib.pyplot as plt
from Bio import SeqIO

GENETIC_CODE = {
    "UUU":"F","UUC":"F","UUA":"L","UUG":"L",
    "UCU":"S","UCC":"S","UCA":"S","UCG":"S",
    "UAU":"Y","UAC":"Y","UAA":"","UAG":"",
    "UGU":"C","UGC":"C","UGA":"*","UGG":"W",
    "CUU":"L","CUC":"L","CUA":"L","CUG":"L",
    "CCU":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAU":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGU":"R","CGC":"R","CGA":"R","CGG":"R",
    "AUU":"I","AUC":"I","AUA":"I","AUG":"M",
    "ACU":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAU":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGU":"S","AGC":"S","AGA":"R","AGG":"R",
    "GUU":"V","GUC":"V","GUA":"V","GUG":"V",
    "GCU":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAU":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGU":"G","GGC":"G","GGA":"G","GGG":"G"
}

AA3 = {
    "A":"Ala","R":"Arg","N":"Asn","D":"Asp","C":"Cys",
    "Q":"Gln","E":"Glu","G":"Gly","H":"His","I":"Ile",
    "L":"Leu","K":"Lys","M":"Met","F":"Phe","P":"Pro",
    "S":"Ser","T":"Thr","W":"Trp","Y":"Tyr","V":"Val",
    "*":"Stop"
}

def clean_dna(seq):
    return "".join(re.findall(r"[ACGT]", seq.upper()))

def dna_to_rna(dna):
    return dna.replace("T", "U")

def count_codons(dna):
    rna = dna_to_rna(dna)
    codons = [rna[i:i+3] for i in range(0, len(rna)-2, 3)]
    return Counter(codons)

def translate(dna):
    rna = dna_to_rna(dna)
    aa = [GENETIC_CODE.get(rna[i:i+3], "?") for i in range(0, len(rna)-2, 3)]
    aa = [a for a in aa if a not in ["*", "?"]]
    return aa

def plot_top_codons(freqs, title):
    top10 = freqs.most_common(10)
    codons, counts = zip(*top10)
    plt.figure(figsize=(8,5))
    plt.bar(codons, counts, color="teal")
    plt.title(title)
    plt.xlabel("Codon")
    plt.ylabel("Frequency")
    plt.show()

def analyze_genome(fasta_path, name):
    record = next(SeqIO.parse(fasta_path, "fasta"))
    dna = clean_dna(str(record.seq))
    codon_freqs = count_codons(dna)
    aa_seq = translate(dna)
    aa_freqs = Counter(aa_seq)
    top3aa = [AA3[a] for a, _ in aa_freqs.most_common(3)]
    print(f"\nTop 3 aminoacids in {name}: {top3aa}")
    plot_top_codons(codon_freqs, f"Top 10 codons in {name}")
    return codon_freqs, aa_freqs, top3aa

# Main execution
covid_codons, covid_aa, covid_top3 = analyze_genome("covid19_genome.fasta", "COVID-19")
flu_codons, flu_aa, flu_top3 = analyze_genome("influenza_genome.fasta", "Influenza")

# Compare common frequent codons
covid_top = set([c for c, _ in covid_codons.most_common(10)])
flu_top = set([c for c, _ in flu_codons.most_common(10)])
common = covid_top & flu_top
print("\nCommon frequent codons between COVID-19 and Influenza:", common)

# AI prompt suggestion
prompt = f"""
Suggest natural foods that are rich in these amino acids:
COVID-19 dominant: {', '.join(covid_top3)}
Influenza dominant: {', '.join(flu_top3)}
List separate plant and animal sources for each.
"""
print("\nAI Prompt:\n", prompt)