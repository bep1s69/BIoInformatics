import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np

class BioinformaticsApp:
    def __init__(self, root):
        self.root = root
        self.root.title("DNA Sequence Analyzer - Sliding Window Analysis")
        self.root.geometry("900x700")
        
        # Variables
        self.sequence = ""
        self.results_data = []
        
        # gui
        self.create_widgets()
        
    def create_widgets(self):
        # Title label
        title_label = tk.Label(self.root, text="DNA Sequence Analyzer", 
                              font=("Arial", 16, "bold"))
        title_label.pack(pady=10)
        
        file_frame = tk.Frame(self.root)
        file_frame.pack(pady=10)
        
        tk.Label(file_frame, text="Select FASTA file:").pack(side=tk.LEFT)
        self.file_path_var = tk.StringVar(value="No file selected")
        tk.Label(file_frame, textvariable=self.file_path_var, 
                fg="blue").pack(side=tk.LEFT, padx=10)
        
        tk.Button(file_frame, text="Browse", 
                 command=self.select_file).pack(side=tk.LEFT, padx=5)
        
        tk.Button(self.root, text="Analyze Sequence", 
                 command=self.analyze_sequence,
                 bg="green", fg="white", font=("Arial", 12)).pack(pady=10)
        
        tk.Button(self.root, text="Export Results to CSV", 
                 command=self.export_results,
                 bg="blue", fg="white", font=("Arial", 10)).pack(pady=5)
        
        # Results display frame
        self.results_frame = tk.Frame(self.root)
        self.results_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Create notebook for tabs
        self.notebook = ttk.Notebook(self.results_frame)
        
        # Tab 1: Data table
        self.table_frame = tk.Frame(self.notebook)
        self.notebook.add(self.table_frame, text="Data Table")
        
        # Tab 2: Chart
        self.chart_frame = tk.Frame(self.notebook)
        self.notebook.add(self.chart_frame, text="Frequency Chart")
        
        self.notebook.pack(fill=tk.BOTH, expand=True)
        
    def select_file(self):
        """Open file dialog to select FASTA file"""
        file_path = filedialog.askopenfilename(
            title="Select FASTA file",
            filetypes=[("FASTA files", "*.fasta *.fa *.fas"), 
                      ("Text files", "*.txt"), 
                      ("All files", "*.*")]
        )
        
        if file_path:
            self.file_path_var.set(file_path.split("/")[-1])  # Show filename only
            self.sequence = self.read_fasta_file(file_path)
            
    def read_fasta_file(self, file_path):
        """Read and work on FASTA file"""
        try:
            sequence = ""
            with open(file_path, 'r') as file:
                for line in file:
                    line = line.strip()
                    if not line.startswith('>'):  
                        sequence += line.upper()  
            
            if not sequence:
                messagebox.showerror("Error", "No sequence found in file!")
                return ""
                
            messagebox.showinfo("Success", f"Loaded sequence with {len(sequence)} bases")
            return sequence
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to read file: {str(e)}")
            return ""
    
    def sliding_window_analysis(self, sequence, window_size=15):
        """Perform sliding window analysis on DNA sequence"""
        if len(sequence) < window_size:
            messagebox.showerror("Error", f"Sequence too short! Need at least {window_size} bases.")
            return []
        
        results = []
        nucleotides = ['A', 'T', 'C', 'G', 'N']  # Include N for ambiguous bases
        
        for i in range(len(sequence) - window_size + 1):
            window = sequence[i:i + window_size]
            
            counts = {nuc: window.count(nuc) for nuc in nucleotides}
            
            frequencies = {nuc: count / window_size for nuc, count in counts.items()}
            
            #store the results
            result = {
                'Window_Start': i + 1,
                'Window_End': i + window_size,
                'A_freq': frequencies['A'],
                'T_freq': frequencies['T'],
                'C_freq': frequencies['C'],
                'G_freq': frequencies['G'],
                'N_freq': frequencies['N'],
                'GC_content': frequencies['G'] + frequencies['C']
            }
            results.append(result)
            
        return results
    
    def analyze_sequence(self):
        """Main analysis function"""
        if not self.sequence:
            messagebox.showerror("Error", "Please select a FASTA file first!")
            return
        
        # Perform sliding window analysis
        self.results_data = self.sliding_window_analysis(self.sequence)
        
        if not self.results_data:
            return
        
        # Display results
        self.display_table()
        self.display_chart()
        
        messagebox.showinfo("Analysis Complete", 
                           f"Analysis completed! {len(self.results_data)} windows analyzed.")
        
    def display_table(self):
        """Display results in table format"""
        # clear previous table
        for widget in self.table_frame.winfo_children():
            widget.destroy()
            
        columns = ['Window_Start', 'Window_End', 'A_freq', 'T_freq', 'C_freq', 'G_freq', 'GC_content']
        
        tree = ttk.Treeview(self.table_frame, columns=columns, show='headings', height=15)
        
        #define headings
        tree.heading('Window_Start', text='Start')
        tree.heading('Window_End', text='End')
        tree.heading('A_freq', text='A Freq')
        tree.heading('T_freq', text='T Freq')
        tree.heading('C_freq', text='C Freq')
        tree.heading('G_freq', text='G Freq')
        tree.heading('GC_content', text='GC Content')
        
        for col in columns:
            tree.column(col, width=90, anchor='center')
        
        display_count = min(100, len(self.results_data))
        for i in range(display_count):
            row = self.results_data[i]
            tree.insert('', 'end', values=[
                row['Window_Start'],
                row['Window_End'],
                f"{row['A_freq']:.3f}",
                f"{row['T_freq']:.3f}",
                f"{row['C_freq']:.3f}",
                f"{row['G_freq']:.3f}",
                f"{row['GC_content']:.3f}"
            ])
        
        v_scrollbar = ttk.Scrollbar(self.table_frame, orient='vertical', command=tree.yview)
        h_scrollbar = ttk.Scrollbar(self.table_frame, orient='horizontal', command=tree.xview)
        tree.configure(yscrollcommand=v_scrollbar.set, xscrollcommand=h_scrollbar.set)
        
        tree.grid(row=0, column=0, sticky='nsew')
        v_scrollbar.grid(row=0, column=1, sticky='ns')
        h_scrollbar.grid(row=1, column=0, sticky='ew')
        
        self.table_frame.grid_rowconfigure(0, weight=1)
        self.table_frame.grid_columnconfigure(0, weight=1)
        
        # Add info label
        info_label = tk.Label(self.table_frame, 
                             text=f"Showing first {display_count} windows of {len(self.results_data)} total",
                             font=("Arial", 10))
        info_label.grid(row=2, column=0, pady=5)
    
    def display_chart(self):
        """Display frequency chart"""
        for widget in self.chart_frame.winfo_children():
            widget.destroy()
            
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
        
        # extract data for plotting
        window_positions = [row['Window_Start'] for row in self.results_data]
        a_freqs = [row['A_freq'] for row in self.results_data]
        t_freqs = [row['T_freq'] for row in self.results_data]
        c_freqs = [row['C_freq'] for row in self.results_data]
        g_freqs = [row['G_freq'] for row in self.results_data]
        gc_content = [row['GC_content'] for row in self.results_data]
        
        # Plot 1: Individual nucleotide frequencies
        ax1.plot(window_positions, a_freqs, label='A', color='red', linewidth=1, alpha=0.8)
        ax1.plot(window_positions, t_freqs, label='T', color='blue', linewidth=1, alpha=0.8)
        ax1.plot(window_positions, c_freqs, label='C', color='green', linewidth=1, alpha=0.8)
        ax1.plot(window_positions, g_freqs, label='G', color='orange', linewidth=1, alpha=0.8)
        
        ax1.set_xlabel('Window Start Position')
        ax1.set_ylabel('Relative Frequency')
        ax1.set_title('DNA Nucleotide Frequencies - Sliding Window Analysis (Window Size: 15)')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim(0, 1)
        
        # Plot 2: GC content
        ax2.plot(window_positions, gc_content, label='GC Content', color='purple', linewidth=1.5)
        ax2.axhline(y=0.5, color='gray', linestyle='--', alpha=0.7, label='50% GC')
        
        ax2.set_xlabel('Window Start Position')
        ax2.set_ylabel('GC Content')
        ax2.set_title('GC Content Along Sequence')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_ylim(0, 1)
        
        plt.tight_layout()
        
        # Embed plot in tkinter
        canvas = FigureCanvasTkAgg(fig, self.chart_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
    def export_results(self):
        """Export results to CSV file"""
        if not self.results_data:
            messagebox.showerror("Error", "No data to export! Please analyze a sequence first.")
            return
            
        file_path = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
            title="Save results as CSV"
        )
        
        if file_path:
            try:
                with open(file_path, 'w') as f:
                    # Write header
                    f.write("Window_Start,Window_End,A_freq,T_freq,C_freq,G_freq,N_freq,GC_content\n")
                    
                    # Write data
                    for row in self.results_data:
                        f.write(f"{row['Window_Start']},{row['Window_End']},"
                               f"{row['A_freq']:.6f},{row['T_freq']:.6f},"
                               f"{row['C_freq']:.6f},{row['G_freq']:.6f},"
                               f"{row['N_freq']:.6f},{row['GC_content']:.6f}\n")
                
                messagebox.showinfo("Export Successful", f"Results exported to {file_path}")
                
            except Exception as e:
                messagebox.showerror("Export Error", f"Failed to export results: {str(e)}")

def main():
    """ main func """
    root = tk.Tk()
    app = BioinformaticsApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()