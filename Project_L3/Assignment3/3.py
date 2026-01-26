import matplotlib.pyplot as plt

 
file_name = 'sequence.fasta'
sequence = ''
window_size = 9

try:
    with open(file_name, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                sequence += line.strip()
except FileNotFoundError:
    print(f"Error: The file '{file_name}' was not found. Please make sure it's in the same directory.")
    exit()

if not sequence:
    print("Error: The FASTA file is empty or does not contain a valid sequence.")
    exit()

positions = []
basic_tm_list = []
salt_tm_list = []

for i in range(len(sequence) - window_size + 1):
    window = sequence[i:i+window_size]
    
    a_count = window.count('A')
    t_count = window.count('T')
    g_count = window.count('G')
    c_count = window.count('C')
    
    # Formula 1: Basic (2+4 rule)
    basic_tm = (a_count + t_count) * 2 + (g_count + c_count) * 4
    
    # Formula 2: Salt-Adjusted
    gc_content = ((g_count + c_count) / window_size) * 100
    salt_tm = 81.5 + (0.41 * gc_content) - (500 / window_size)
    
    positions.append(i)
    basic_tm_list.append(basic_tm)
    salt_tm_list.append(salt_tm)

 
print("--- Analysis Complete ---")
print(f"Basic Formula (2+4 Rule) - Min Tm: {min(basic_tm_list):.2f}°C, Max Tm: {max(basic_tm_list):.2f}°C")
print(f"Salt-Adjusted Formula    - Min Tm: {min(salt_tm_list):.2f}°C, Max Tm: {max(salt_tm_list):.2f}°C")
print("-" * 25)

# --- 3. NEW: Get User Input for Threshold ---
threshold = None
while threshold is None:
    try:
        threshold_input = input("Enter the melting temperature threshold (°C) to filter results: ")
        threshold = float(threshold_input)
    except ValueError:
        print("Invalid input. Please enter a number (e.g., 35.5).")

 
basic_bars = [(i, 1) for i, tm in enumerate(basic_tm_list) if tm > threshold]
salt_bars = [(i, 1) for i, tm in enumerate(salt_tm_list) if tm > threshold]
 
 
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 10), sharex=True)
fig.suptitle('DNA Melting Temperature Analysis', fontsize=16)

 
ax1.plot(positions, basic_tm_list, label='Basic Formula (2+4 Rule)', color='blue')
ax1.plot(positions, salt_tm_list, label='Salt-Adjusted Formula', color='orangered')
ax1.axhline(y=threshold, color='green', linestyle='--', label=f'Threshold ({threshold}°C)')
ax1.set_title('Melting Temperature Along Sequence')
ax1.set_ylabel('Melting Temperature (Tm) in °C')
ax1.grid(True)
ax1.legend()

 
ax2.set_title(f'Positions with Tm > {threshold}°C')

 
 
ax2.broken_barh(basic_bars, (10, 8), facecolors='blue')
ax2.broken_barh(salt_bars, (22, 8), facecolors='green')

 
ax2.set_xlabel('Starting Position of 9-base Window')
ax2.set_yticks([14, 26])  
ax2.set_yticklabels(['Basic Formula', 'Salt-Adjusted Formula'])
ax2.set_ylim(5, 35) 
ax2.grid(axis='x')  

 
plt.tight_layout(rect=[0, 0, 1, 0.96])  
plt.show()