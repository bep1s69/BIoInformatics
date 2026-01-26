import os
from Bio import SeqIO
from Bio.Seq import Seq

# --- CONFIGURATION BASED ON ASSIGNMENT ---
# Constraints from the blackboard
MIN_IR_LENGTH = 4
MAX_IR_LENGTH = 6
# Maximum allowed length for the entire Transposable Element (TE) body
MAX_TE_LENGTH = 5000  # Example max size in base pairs (5 kb)

# --- HELPER FUNCTIONS ---

def find_inverted_repeats(sequence_record):
    """
    Scans a sequence for pairs of Inverted Repeats (IRs) within a defined distance.
    """
    results = []
    seq = str(sequence_record.seq)
    genome_length = len(seq)
    
    print(f"-> Analyzing sequence: {sequence_record.id} (Length: {genome_length})")
    
    # 1. Iterate through all possible start positions (i) for the first IR (IR1)
    for i in range(genome_length):
        
        # 2. Iterate through all possible lengths for IR1
        for ir_len in range(MIN_IR_LENGTH, MAX_IR_LENGTH + 1):
            
            # Check if IR1 fits within the genome bounds
            if i + ir_len > genome_length:
                continue
                
            ir1_sequence = seq[i : i + ir_len]
            # Calculate the expected sequence for the second IR (IR2)
            # IR2 must be the Reverse Complement of IR1
            expected_ir2 = str(Seq(ir1_sequence).reverse_complement())
            
            # 3. Search for IR2 within the valid distance
            # The search starts immediately after IR1 ends (i + ir_len)
            # and must end before the maximum TE length is exceeded.
            start_search = i + ir_len
            end_search = min(genome_length, i + ir_len + MAX_TE_LENGTH)
            
            # Iterate through all possible start positions (j) for the second IR (IR2)
            for j in range(start_search, end_search):
                
                # Check if IR2 fits within the genome bounds
                if j + ir_len > genome_length:
                    continue
                    
                ir2_sequence = seq[j : j + ir_len]
                
                # 4. Check for a match
                if ir2_sequence == expected_ir2:
                    # Found a potential Transposable Element!
                    
                    # Position of TE: Start of IR1 to End of IR2
                    te_start = i + 1  # 1-based indexing for output
                    te_end = j + ir_len
                    te_length = te_end - te_start + 1
                    
                    # The TE element body is the sequence between the two IRs,
                    # but the output asks for the position and length of the element # which typically includes the IRs.
                    
                    results.append({
                        'genome_id': sequence_record.id,
                        'position_start': te_start,
                        'position_end': te_end,
                        'length': te_length,
                        'ir_length': ir_len,
                        'ir1_seq': ir1_sequence
                    })
                    
    return results

def process_genomes(genome_files):
    """
    Loads and processes a list of genome FASTA files.
    """
    all_results = []
    
    for file_path in genome_files:
        print(f"\n--- Processing File: {file_path} ---")
        
        # Biopython's SeqIO handles stream processing efficiently
        # We assume one record (chromosome) per file for typical bacterial genomes
        try:
            for record in SeqIO.parse(file_path, "fasta"):
                te_detections = find_inverted_repeats(record)
                all_results.extend(te_detections)
                
        except FileNotFoundError:
            print(f"ERROR: File not found at {file_path}. Skipping.")
        except Exception as e:
            print(f"An error occurred while processing {file_path}: {e}")
            
    return all_results

def display_results(results):
    """
    Displays the detected TEs in a clear, formatted table.
    """
    if not results:
        print("\n--- RESULTS ---")
        print("No potential Transposable Elements detected based on criteria.")
        return

    print("\n--- DETECTED TRANSPOSEABLE ELEMENTS (Based on IR Structure) ---")
    print(f"{'Genome ID':<20} | {'Start':<8} | {'End':<8} | {'Length (bp)':<12} | {'IR Length':<10} | {'IR1 Sequence':<15}")
    print("-" * 80)
    for res in results:
        print(
            f"{res['genome_id'][:20]:<20} | {res['position_start']:<8} | {res['position_end']:<8} | {res['length']:<12} | {res['ir_length']:<10} | {res['ir1_seq']:<15}"
        )


# --- MAIN EXECUTION ---
if __name__ == "__main__":
    
    # ⚠ 1. IMPORTANT: Replace these paths with the actual locations 
    # of the three bacterial FASTA files you downloaded from NCBI.
    # For testing, you must create a placeholder file first.
    
    # Example placeholder paths - CHANGE THIS:
    GENOME_PATHS = [
        "path/to/your/genome1.fasta",
        "path/to/your/genome2.fasta",
        "path/to/your/genome3.fasta",
    ]

    # --- SIMULATION (REMOVE AFTER REAL FILES ARE DOWNLOADED) ---
    # Create a small test file for demonstration if the real files aren't ready
    try:
        os.makedirs("test_data", exist_ok=True)
        TEST_FILE = "test_data/test_genome.fasta"
        # Example sequence with a 5bp IR (TTAAG) and a 4bp IR (GTCA)
        # Element 1 (Length: 10 + 5*2 = 20): TTAAG...C T T T T T T T T T C T T A A
        #                       Reverse complement of TTAAG is CTTAA
        # Element 2 (Length: 10 + 4*2 = 18): GTCA...G T T G T G T G T T T G A C
        #                       Reverse complement of GTCA is TGAC
        test_seq = ">E_coli_K12\n" \
                   "ATATAC*TTAAGCTTTTTTTTTTCTTAAGATGATGGTCAGTTGTGTGTGTTGAC*AAAAA"
        with open(TEST_FILE, "w") as f:
            f.write(test_seq)
        
        GENOME_PATHS = [TEST_FILE] # Use only the test file for initial run
        
    except Exception as e:
        print(f"Could not create test file: {e}")
    # --- END OF SIMULATION ---

    # Process the files and get the results
    detected_elements = process_genomes(GENOME_PATHS)
    
    # Display the final output
    display_results(detected_elements)