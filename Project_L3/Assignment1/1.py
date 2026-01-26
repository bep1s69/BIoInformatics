import math
import sys

def calculate_tm_simple(sequence):
    sequence = sequence.upper().replace(' ', '').replace('\n', '').replace('\r', '')
    if not sequence or not all(c in 'ATGC' for c in sequence):
        return None
    g_c_count = sequence.count('G') + sequence.count('C')
    return 4 * g_c_count + 2 * (len(sequence) - g_c_count)

def calculate_tm_advanced(sequence, na_conc=0.05):
    sequence = sequence.upper().replace(' ', '').replace('\n', '').replace('\r', '')
    if not sequence or not all(c in 'ATGC' for c in sequence):
        return None
    length = len(sequence)
    gc_percent = 100 * (sequence.count('G') + sequence.count('C')) / length
    log_na = math.log10(na_conc)
    tm = 81.5 + 16.6 * log_na + 0.41 * gc_percent - 600 / length
    return tm

def main():
    print("DNA Melting Temperature Calculator")
     
    print("Enter a DNA sequence using only A, T, G, or C ")
    print("The program will calculate the melting temperature using two formulas.")
    
    try:
        sequence = input("Enter DNA sequence: ").strip()
    except EOFError:
        print("Error: Input was interrupted. Please provide a valid DNA sequence.", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("Error: Input was cancelled. Program exiting.", file=sys.stderr)
        sys.exit(1)
    print("\nInput Sequence:", sequence)
    
    tm_simple = calculate_tm_simple(sequence)
    if tm_simple is not None:
        print(f"Simple Formula: {tm_simple} °C")
        print("  Formula: Tm = 4 * (G + C) + 2 * (A + T)")
    else:
        print("Error: Invalid DNA sequence. Use only A, T, G, C.", file=sys.stderr)
        sys.exit(1)
    tm_adv = calculate_tm_advanced(sequence)
    if tm_adv is not None:
        print(f"Advanced Formula: {tm_adv:.1f} °C")
        print("  Formula: Tm = 81.5 + 16.6 * log10([Na+]) + 0.41 * (%GC) - 600 / length")
        print("  Using [Na+] = 50 mM")
    else:
        print("Error: Invalid DNA sequence. Use only A, T, G, C.", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()