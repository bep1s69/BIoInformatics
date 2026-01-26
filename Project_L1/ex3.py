import React, { useState } from 'react';
import { Upload, FileText, BarChart3, Search } from 'lucide-react';

const FastaAnalyzer = () => {
  const [fastaContent, setFastaContent] = useState('');
  const [results, setResults] = useState('');
  const [fileName, setFileName] = useState('');

  // Algorithm 1: Calculate sequence statistics
  const calculateStats = (sequence) => {
    const seq = sequence.toUpperCase();
    const length = seq.length;
    
    const counts = {
      A: (seq.match(/A/g) || []).length,
      T: (seq.match(/T/g) || []).length,
      G: (seq.match(/G/g) || []).length,
      C: (seq.match(/C/g) || []).length,
      U: (seq.match(/U/g) || []).length
    };
    
    const gcContent = length > 0 ? ((counts.G + counts.C) / length * 100).toFixed(2) : 0;
    
    return {
      length,
      counts,
      gcContent,
      atContent: (100 - gcContent).toFixed(2)
    };
  };

  // Algorithm 2: Find repeating patterns (motifs)
  const findMotifs = (sequence, minLength = 3, maxLength = 6) => {
    const seq = sequence.toUpperCase();
    const motifs = {};
    
    for (let len = minLength; len <= maxLength; len++) {
      for (let i = 0; i <= seq.length - len; i++) {
        const motif = seq.substring(i, i + len);
        if (!/[^ATGCU]/.test(motif)) {
          motifs[motif] = (motifs[motif] || 0) + 1;
        }
      }
    }
    
    return Object.entries(motifs)
      .filter(([_, count]) => count > 1)
      .sort((a, b) => b[1] - a[1])
      .slice(0, 10);
  };

  // Parse FASTA file
  const parseFasta = (content) => {
    const lines = content.split('\n');
    let header = '';
    let sequence = '';
    
    for (let line of lines) {
      line = line.trim();
      if (line.startsWith('>')) {
        header = line.substring(1);
      } else if (line.length > 0) {
        sequence += line;
      }
    }
    
    return { header, sequence };
  };

  // Process FASTA file
  const processFasta = () => {
    if (!fastaContent) {
      setResults('Please load a FASTA file first.');
      return;
    }

    const { header, sequence } = parseFasta(fastaContent);
    
    if (!sequence) {
      setResults('Error: No valid sequence found in file.');
      return;
    }

    // Run Algorithm 1: Statistics
    const stats = calculateStats(sequence);
    
    // Run Algorithm 2: Motif finding
    const motifs = findMotifs(sequence);
    
    // Format results
    let output = '═══════════════════════════════════════════════\n';
    output += '          FASTA FILE ANALYSIS RESULTS\n';
    output += '═══════════════════════════════════════════════\n\n';
    
    output += ' SEQUENCE INFORMATION:\n';
    output += '─────────────────────────────────────────────\n';
    output += Header: ${header}\n;
    output += Sequence Length: ${stats.length} bp\n\n;
    
    output += ' ALGORITHM 1 - NUCLEOTIDE COMPOSITION:\n';
    output += '─────────────────────────────────────────────\n';
    output += Adenine (A):  ${stats.counts.A} (${(stats.counts.A/stats.length*100).toFixed(2)}%)\n;
    output += Thymine (T):  ${stats.counts.T} (${(stats.counts.T/stats.length*100).toFixed(2)}%)\n;
    output += Guanine (G):  ${stats.counts.G} (${(stats.counts.G/stats.length*100).toFixed(2)}%)\n;
    output += Cytosine (C): ${stats.counts.C} (${(stats.counts.C/stats.length*100).toFixed(2)}%)\n;
    output += Uracil (U):   ${stats.counts.U} (${(stats.counts.U/stats.length*100).toFixed(2)}%)\n\n;
    
    output += GC Content: ${stats.gcContent}%\n;
    output += AT Content: ${stats.atContent}%\n\n;
    
    output += ' ALGORITHM 2 - REPEATING MOTIFS:\n';
    output += '─────────────────────────────────────────────\n';
    if (motifs.length > 0) {
      output += 'Top 10 repeating patterns:\n\n';
      motifs.forEach(([motif, count], idx) => {
        output += ${idx + 1}. "${motif}" - occurs ${count} times\n;
      });
    } else {
      output += 'No significant repeating patterns found.\n';
    }
    
    output += '\n═══════════════════════════════════════════════\n';
    
    setResults(output);
  };

  // Handle file upload
  const handleFileUpload = async (e) => {
    const file = e.target.files[0];
    if (file) {
      setFileName(file.name);
      const content = await window.fs.readFile(file.name, { encoding: 'utf8' });
      setFastaContent(content);
      setResults('File loaded successfully! Click "Analyze Sequence" to process.');
    }
  };

  // Load sample FASTA
  const loadSampleFasta = () => {
    const sample = `>SEQ001 | Custom synthetic DNA sequence | Length: 12bp | GC-rich region
ATTGCCCCGAAT`;
    
    setFastaContent(sample);
    setFileName('custom_sequence.fasta');
    setResults('Custom DNA sequence "ATTGCCCCGAAT" loaded! Click "Analyze Sequence" to process.');
  };

  return (
    <div className="min-h-screen bg-gradient-to-br from-blue-50 to-indigo-100 p-8">
      <div className="max-w-6xl mx-auto">
        <div className="bg-white rounded-lg shadow-xl p-8">
          <div className="flex items-center gap-3 mb-6">
            <FileText className="w-8 h-8 text-indigo-600" />
            <h1 className="text-3xl font-bold text-gray-800">FASTA Sequence Analyzer</h1>
          </div>
          
          <div className="mb-6 p-4 bg-blue-50 rounded-lg border border-blue-200">
            <p className="text-sm text-gray-700">
              <strong>Features:</strong> Upload or load a sample FASTA file to analyze DNA/RNA sequences.
              The application runs two algorithms: (1) Nucleotide composition analysis, and (2) Repeating motif detection.
            </p>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-4 mb-6">
            <div className="flex flex-col gap-3">
              <label className="flex items-center justify-center gap-2 px-6 py-3 bg-indigo-600 text-white rounded-lg hover:bg-indigo-700 cursor-pointer transition">
                <Upload className="w-5 h-5" />
                <span>Choose FASTA File</span>
                <input 
                  type="file" 
                  accept=".fasta,.fa,.fna,.txt"
                  onChange={handleFileUpload}
                  className="hidden"
                />
              </label>
              
              <button
                onClick={loadSampleFasta}
                className="flex items-center justify-center gap-2 px-6 py-3 bg-green-600 text-white rounded-lg hover:bg-green-700 transition"
              >
                <FileText className="w-5 h-5" />
                <span>Load Sample FASTA</span>
              </button>
            </div>

            <button
              onClick={processFasta}
              className="flex items-center justify-center gap-2 px-6 py-3 bg-purple-600 text-white rounded-lg hover:bg-purple-700 transition disabled:bg-gray-400 disabled:cursor-not-allowed"
              disabled={!fastaContent}
            >
              <BarChart3 className="w-5 h-5" />
              <span>Analyze Sequence</span>
            </button>
          </div>

          {fileName && (
            <div className="mb-4 p-3 bg-gray-100 rounded-lg">
              <p className="text-sm text-gray-700">
                <strong>Loaded file:</strong> {fileName}
              </p>
            </div>
          )}

          <div className="bg-gray-50 rounded-lg p-6 border-2 border-gray-200">
            <h2 className="text-lg font-semibold text-gray-800 mb-3 flex items-center gap-2">
              <Search className="w-5 h-5" />
              Analysis Results
            </h2>
            <textarea
              value={results}
              readOnly
              placeholder="Results will appear here after analysis..."
              className="w-full h-96 p-4 bg-white border border-gray-300 rounded-lg font-mono text-sm resize-none focus:outline-none focus:ring-2 focus:ring-indigo-500"
            />
          </div>
        </div>
      </div>
    </div>
  );
};

export default FastaAnalyzer;