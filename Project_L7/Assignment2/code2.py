import time
import matplotlib.pyplot as plt
from Bio import Entrez, SeqIO
from collections import defaultdict

 
def find_tandem_repeats(sequence, min_len=3, max_len=6, min_repeats=2):
    """
    Finds all tandem repeats in a DNA sequence for motifs of a given length.
    """
    
    # Clean the sequence: remove newlines/spaces and convert to uppercase
    sequence = "".join(str(sequence).split()).upper()
    seq_len = len(sequence)
    results = {}

    for k in range(min_len, max_len + 1):
        i = 0
        while i <= seq_len - k:
            motif = sequence[i : i + k]
            
            # Skip motifs with non-standard DNA/RNA characters (like 'N')
            if not all(c in 'ATCGU' for c in motif):
                i += 1
                continue

            repeat_count = 1
            
            while True:
                next_start = i + (repeat_count * k)
                if next_start + k > seq_len:
                    break
                
                next_motif = sequence[next_start : next_start + k]
                
                if next_motif == motif:
                    repeat_count += 1
                else:
                    break
            
            if repeat_count >= min_repeats:
                if motif not in results:
                    results[motif] = []
                    
                result_entry = {
                    'position': i,
                    'count': repeat_count,
                    'full_sequence': sequence[i : i + (repeat_count * k)]
                }
                results[motif].append(result_entry)
                
                # Skip past the entire repeat block we just found
                i += (repeat_count * k)
            else:
                # No repeat, just move to the next base
                i += 1
                
    return results

# --- PART 2: DOWNLOAD AND ANALYSIS ---

def download_and_analyze_genomes(accession_ids):
    """
    Downloads, analyzes, and summarizes repeat data for a list of NCBI accessions.
    """
    
    # You must provide your email to use NCBI's API
    Entrez.email = "your_email@example.com" 
    
    # This list will store our final summary data for plotting
    plot_data = []

    print("--- Starting Genome Analysis ---")
    
    for accession_id in accession_ids:
        print(f"\n[1/3] Downloading: {accession_id}...")
        
        try:
            # Fetch the sequence from NCBI Nucleotide database in FASTA format
            handle = Entrez.efetch(db="nucleotide", 
                                   id=accession_id, 
                                   rettype="fasta", 
                                   retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()
            
            sequence = record.seq
            print(f"[2/3] Analyzing {accession_id} (Length: {len(sequence)} bp)...")

            # Run the repeat-finding algorithm
            repeats = find_tandem_repeats(sequence)
            
            if not repeats:
                print(f"[3/3] No repeats found for {accession_id}.")
                plot_data.append({'id': accession_id, 'motif': 'None', 'frequency': 0})
                continue

            # --- Find the MOST frequent motif ---
            # "Most frequent" is defined as the motif that appears
            # in the most locations (not the longest single repeat).
            
            most_frequent_motif = ""
            max_frequency = 0
            
            for motif, locations in repeats.items():
                current_frequency = len(locations)
                if current_frequency > max_frequency:
                    max_frequency = current_frequency
                    most_frequent_motif = motif
            
            print(f"[3/3] Complete. Top repeat: '{most_frequent_motif}' (found at {max_frequency} locations).")
            
            # Save the summary for plotting
            plot_data.append({'id': accession_id, 'motif': most_frequent_motif, 'frequency': max_frequency})
            
            # Be polite to the NCBI server
            time.sleep(1) 

        except Exception as e:
            print(f"  Error processing {accession_id}: {e}")
            
    return plot_data

 

def plot_repeat_frequenies(plot_data):
    """
    Uses Matplotlib to create a bar chart of the most frequent repeats.
    """
    if not plot_data:
        print("No data to plot.")
        return
 
    genome_ids = [entry['id'] for entry in plot_data]
    frequencies = [entry['frequency'] for entry in plot_data]
    motifs = [entry['motif'] for entry in plot_data]

    plt.figure(figsize=(15, 8)) 
    
    bars = plt.bar(genome_ids, frequencies, color='skyblue')
    
    plt.title('Most Frequent Tandem Repeat (3-6bp) in 10 Flu Genomes', fontsize=16)
    plt.ylabel('Number of Repeat Locations', fontsize=12)
    plt.xlabel('Genome Segment (NCBI Accession ID)', fontsize=12)
    
   
    plt.xticks(rotation=45, ha="right")
    
    
    for i, bar in enumerate(bars):
        yval = bar.get_height()
        motif_label = motifs[i]
        plt.text(bar.get_x() + bar.get_width()/2.0, 
                 yval + 0.1,  # Position label just above the bar
                 motif_label, 
                 ha='center', 
                 va='bottom',
                 fontsize=10,
                 fontweight='bold')
    
    plt.tight_layout() # Adjust layout to fit labels
    print("\n--- Analysis complete. Displaying plot. ---")
    plt.show()

 
if __name__ == "__main__":
    # 1. Define the 10 genomes to analyze
    accession_list = [
        "OQ837085.1", "OQ837084.1", "OQ837088.1", "OQ837087.1",
        "PP407268.1", "PP407269.1", "PP407270.1", "PP407271.1",
        "MG692771.1", "MG692772.1"
    ]
     
    summary_data = download_and_analyze_genomes(accession_list)
    
    
    plot_repeat_frequenies(summary_data)