import React, { useState } from 'react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from 'recharts';
import { Upload, Activity } from 'lucide-react';

const DNAMeltingTempAnalyzer = () => {
  const [data, setData] = useState([]);
  const [seqInfo, setSeqInfo] = useState(null);
  const [error, setError] = useState('');

  // Wallace formula: Tm = 2(A+T) + 4(G+C)
  const calculateWallaceTm = (seq) => {
    const counts = {A: 0, T: 0, G: 0, C: 0};
    for (let base of seq.toUpperCase()) {
      if (counts.hasOwnProperty(base)) counts[base]++;
    }
    return 2 * (counts.A + counts.T) + 4 * (counts.G + counts.C);
  };

  // Nearest-neighbor approximation (simplified)
  // Tm = 64.9 + 41 * (G+C-16.4) / length
  const calculateNearestNeighborTm = (seq) => {
    const upper = seq.toUpperCase();
    let gc = 0;
    for (let base of upper) {
      if (base === 'G' || base === 'C') gc++;
    }
    const gcContent = gc / seq.length;
    return 64.9 + 41 * (gcContent - 0.164);
  };

  const parseFasta = (text) => {
    const lines = text.split('\n').filter(line => line.trim());
    let header = '';
    let sequence = '';
    
    for (let line of lines) {
      if (line.startsWith('>')) {
        header = line.substring(1).trim();
      } else {
        sequence += line.trim().replace(/\s/g, '');
      }
    }
    
    return { header, sequence };
  };

  const analyzeDNA = (sequence) => {
    const windowSize = 9;
    const results = [];
    
    if (sequence.length < windowSize) {
      setError(`Sequence too short. Need at least ${windowSize} bases.`);
      return [];
    }

    for (let i = 0; i <= sequence.length - windowSize; i++) {
      const window = sequence.substring(i, i + windowSize);
      const wallaceTm = calculateWallaceTm(window);
      const nnTm = calculateNearestNeighborTm(window);
      
      results.push({
        position: i + 1,
        window: window,
        wallace: parseFloat(wallaceTm.toFixed(2)),
        nearestNeighbor: parseFloat(nnTm.toFixed(2))
      });
    }
    
    return results;
  };

  const handleFileUpload = async (event) => {
    const file = event.target.files[0];
    if (!file) return;

    setError('');
    setData([]);
    setSeqInfo(null);

    try {
      const text = await file.text();
      const { header, sequence } = parseFasta(text);
      
      if (!sequence) {
        setError('No valid DNA sequence found in file');
        return;
      }

      // Validate DNA sequence
      const validBases = /^[ATGCN]+$/i;
      if (!validBases.test(sequence)) {
        setError('Invalid DNA sequence. Only A, T, G, C, and N are allowed.');
        return;
      }

      const results = analyzeDNA(sequence);
      
      if (results.length > 0) {
        setData(results);
        setSeqInfo({
          header: header || 'Unnamed sequence',
          length: sequence.length,
          windows: results.length
        });
      }
    } catch (err) {
      setError(`Error reading file: ${err.message}`);
    }
  };

  const handlePasteSequence = () => {
    const seq = prompt('Paste your DNA sequence (FASTA format or raw sequence):');
    if (!seq) return;

    setError('');
    setData([]);
    setSeqInfo(null);

    const { header, sequence } = parseFasta(seq);
    
    if (!sequence) {
      setError('No valid DNA sequence found');
      return;
    }

    const validBases = /^[ATGCN]+$/i;
    if (!validBases.test(sequence)) {
      setError('Invalid DNA sequence. Only A, T, G, C, and N are allowed.');
      return;
    }

    const results = analyzeDNA(sequence);
    
    if (results.length > 0) {
      setData(results);
      setSeqInfo({
        header: header || 'Pasted sequence',
        length: sequence.length,
        windows: results.length
      });
    }
  };

  return (
    <div className="min-h-screen bg-gradient-to-br from-blue-50 to-indigo-100 p-8">
      <div className="max-w-6xl mx-auto">
        <div className="bg-white rounded-lg shadow-xl p-8">
          <div className="flex items-center gap-3 mb-6">
            <Activity className="w-8 h-8 text-indigo-600" />
            <h1 className="text-3xl font-bold text-gray-800">
              DNA Melting Temperature Analyzer
            </h1>
          </div>

          <div className="mb-6 p-4 bg-blue-50 rounded-lg border border-blue-200">
            <h2 className="font-semibold text-blue-900 mb-2">How it works:</h2>
            <ul className="text-sm text-blue-800 space-y-1">
              <li>• Sliding window size: <strong>9 nucleotides</strong></li>
              <li>• <strong>Wallace Formula:</strong> Tm = 2(A+T) + 4(G+C) - Simple, quick approximation</li>
              <li>• <strong>Nearest-Neighbor:</strong> Tm = 64.9 + 41 * (GC% - 0.164) - More accurate for longer oligos</li>
            </ul>
          </div>

          <div className="flex gap-4 mb-6">
            <label className="flex-1 cursor-pointer">
              <div className="flex items-center justify-center gap-2 px-6 py-3 bg-indigo-600 text-white rounded-lg hover:bg-indigo-700 transition-colors">
                <Upload className="w-5 h-5" />
                <span>Upload FASTA File</span>
              </div>
              <input
                type="file"
                accept=".fasta,.fa,.txt"
                onChange={handleFileUpload}
                className="hidden"
              />
            </label>
            
            <button
              onClick={handlePasteSequence}
              className="px-6 py-3 bg-green-600 text-white rounded-lg hover:bg-green-700 transition-colors"
            >
              Paste Sequence
            </button>
          </div>

          {error && (
            <div className="mb-6 p-4 bg-red-50 border border-red-200 rounded-lg text-red-700">
              {error}
            </div>
          )}

          {seqInfo && (
            <div className="mb-6 p-4 bg-gray-50 rounded-lg border border-gray-200">
              <h3 className="font-semibold text-gray-700 mb-2">Sequence Information:</h3>
              <div className="grid grid-cols-3 gap-4 text-sm">
                <div>
                  <span className="text-gray-600">Header:</span>
                  <p className="font-medium text-gray-800 truncate">{seqInfo.header}</p>
                </div>
                <div>
                  <span className="text-gray-600">Length:</span>
                  <p className="font-medium text-gray-800">{seqInfo.length} bp</p>
                </div>
                <div>
                  <span className="text-gray-600">Windows analyzed:</span>
                  <p className="font-medium text-gray-800">{seqInfo.windows}</p>
                </div>
              </div>
            </div>
          )}

          {data.length > 0 && (
            <div className="mt-8">
              <h2 className="text-xl font-semibold text-gray-800 mb-4">
                Melting Temperature Profile
              </h2>
              <ResponsiveContainer width="100%" height={400}>
                <LineChart data={data}>
                  <CartesianGrid strokeDasharray="3 3" />
                  <XAxis 
                    dataKey="position" 
                    label={{ value: 'Position in Sequence', position: 'insideBottom', offset: -5 }}
                  />
                  <YAxis 
                    label={{ value: 'Temperature (°C)', angle: -90, position: 'insideLeft' }}
                  />
                  <Tooltip 
                    content={({ active, payload }) => {
                      if (active && payload && payload.length) {
                        const data = payload[0].payload;
                        return (
                          <div className="bg-white p-3 border border-gray-300 rounded shadow-lg">
                            <p className="font-semibold">Position: {data.position}</p>
                            <p className="text-xs font-mono text-gray-600 mb-2">{data.window}</p>
                            <p className="text-blue-600">Wallace: {data.wallace}°C</p>
                            <p className="text-orange-600">Nearest-Neighbor: {data.nearestNeighbor}°C</p>
                          </div>
                        );
                      }
                      return null;
                    }}
                  />
                  <Legend />
                  <Line 
                    type="monotone" 
                    dataKey="wallace" 
                    stroke="#3b82f6" 
                    name="Wallace Formula"
                    strokeWidth={2}
                    dot={false}
                  />
                  <Line 
                    type="monotone" 
                    dataKey="nearestNeighbor" 
                    stroke="#f97316" 
                    name="Nearest-Neighbor"
                    strokeWidth={2}
                    dot={false}
                  />
                </LineChart>
              </ResponsiveContainer>

              <div className="mt-6 p-4 bg-gray-50 rounded-lg">
                <h3 className="font-semibold text-gray-700 mb-2">Statistics:</h3>
                <div className="grid grid-cols-2 gap-4 text-sm">
                  <div>
                    <p className="text-gray-600">Wallace Tm Range:</p>
                    <p className="font-medium text-blue-600">
                      {Math.min(...data.map(d => d.wallace)).toFixed(2)}°C - {Math.max(...data.map(d => d.wallace)).toFixed(2)}°C
                    </p>
                  </div>
                  <div>
                    <p className="text-gray-600">Nearest-Neighbor Tm Range:</p>
                    <p className="font-medium text-orange-600">
                      {Math.min(...data.map(d => d.nearestNeighbor)).toFixed(2)}°C - {Math.max(...data.map(d => d.nearestNeighbor)).toFixed(2)}°C
                    </p>
                  </div>
                </div>
              </div>
            </div>
          )}

          {data.length === 0 && !error && (
            <div className="text-center py-12 text-gray-500">
              <Activity className="w-16 h-16 mx-auto mb-4 opacity-50" />
              <p>Upload a FASTA file or paste a DNA sequence to begin analysis</p>
            </div>
          )}
        </div>
      </div>
    </div>
  );
};

export default DNAMeltingTempAnalyzer;
