import os
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio import motifs
import matplotlib.pyplot as plt

Entrez.email = "your.mail@example.com" 
GENOME_IDS = [
    "NC_002017", "NC_002018", "NC_002019", "NC_002020", 
    "NC_002021", "NC_002022", "NC_002023", "NC_004908",
    "GQ505938", "GQ505939"
]

#define the motif
instances = [Seq("GATAAA"), Seq("GATACT"), Seq("GATAAA"), Seq("GATAAA")]
m = motifs.create(instances)
pwm = m.counts.normalize(pseudocounts=0.5) # Adding pseudocounts prevents log(0) errors

pssm = pwm.log_odds()

def download_and_scan():
    for gid in GENOME_IDS:
        print(f"Processing {gid}...")
        
        try:
            #download genome
            handle = Entrez.efetch(db="nucleotide", id=gid, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()
            
            sequence = record.seq
            
            #scan for Motifs
            #calculate() returns a Bio.motif.PositionScoreMatrix object
            scores = pssm.calculate(sequence)
            
            plt.figure(figsize=(12, 4))
            plt.plot(scores)
            plt.title(f"Motif Signal: {gid} - {record.description[:40]}")
            plt.xlabel("Position")
            plt.ylabel("Log-odds Score")
            
            #highlight likely functional sites (threshold = 90% of max score)
            threshold = max(scores) * 0.9
            plt.axhline(y=threshold, color='r', linestyle='--', label='Likely Motif')
            
            plt.tight_layout()
            plt.show()
            
        except Exception as e:
            print(f"Error processing {gid}: {e}")

if __name__ == "__main__":
    download_and_scan()