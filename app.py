import os
import io
import time
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
from flask import Flask, request, jsonify, render_template 
from flask_cors import CORS
import argparse  # NEW: For handling permanent CLI arguments
import glob    # NEW: For file path scanning

# --- FLASK SETUP ---
app = Flask(__name__)
CORS(app)

UPLOAD_FOLDER = 'uploads'
STATIC_FOLDER = 'static'

os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(STATIC_FOLDER, exist_ok=True)


@app.route('/')
def home():
    """Renders the HTML interface for the interactive demo."""
    return render_template('i.html')

# --- CORE BIOINFORMATICS FUNCTIONS ---

def load_sequence(fasta_path):
    """Loads and returns the sequence string from a single FASTA file."""
    try:
        with open(fasta_path, "r") as handle:
            record = next(SeqIO.parse(handle, "fasta"))
            return str(record.seq).upper()
    except Exception as e:
        print(f"Error loading {fasta_path}: {e}")
        return ""

def generate_kmers(sequence, k):
    """Generates a list of k-mers (features) from a given sequence."""
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

def compute_similarity(seq1, seq2, k=21):
    """Calculates k-mer based similarity (Forward and Jaccard) between two sequences."""
    kmers1, kmers2 = generate_kmers(seq1, k), generate_kmers(seq2, k)
    if not kmers1 or not kmers2: return 0, 0
    counts1, counts2 = Counter(kmers1), Counter(kmers2)
    common_kmers = set(counts1) & set(counts2)
    
    # Calculate Forward Similarity: Shared kmers / Length of first sequence
    shared = sum(min(counts1[k_val], counts2[k_val]) for k_val in common_kmers)
    forward_sim = (shared / len(kmers1)) * 100 if len(kmers1) > 0 else 0
    
    # Calculate Jaccard Similarity: Intersection / Union
    union_kmers = set(counts1) | set(counts2)
    jaccard_sim = (len(common_kmers) / len(union_kmers)) * 100 if len(union_kmers) > 0 else 0
    return forward_sim, jaccard_sim

def compare_all_genomes(file_paths, k=21):
    """Compares all genomes against each other and returns the result DataFrame."""
    sequences = {}
    for p in file_paths:
        seq = load_sequence(p)
        if seq:
            sequences[os.path.basename(p)] = seq

    filenames = list(sequences.keys())
    results = []
    
    for i in range(len(filenames)):
        for j in range(i + 1, len(filenames)):
            file1, file2 = filenames[i], filenames[j]
            seq1, seq2 = sequences[file1], sequences[file2]
            
            # Use reverse complement for comparison, common in genomics
            seq2_rc = str(Seq(seq2).reverse_complement())
            
            forward_sim, jaccard_sim = compute_similarity(seq1, seq2, k)
            reverse_sim, jaccard_sim_rc = compute_similarity(seq1, seq2_rc, k)
            
            results.append([file1, file2, round(forward_sim, 2), round(reverse_sim, 2), round(jaccard_sim, 2), round(jaccard_sim_rc, 2)])
            
    df = pd.DataFrame(results, columns=["Genome1", "Genome2", "Forward_Similarity(%)", "Reverse_Similarity(%)", "Jaccard_Forward(%)", "Jaccard_Reverse(%)"])
    df = df[["Genome1", "Genome2", "Jaccard_Forward(%)", "Jaccard_Reverse(%)", "Reverse_Similarity(%)", "Forward_Similarity(%)"]]
    return df

def generate_and_save_plots(df, threshold=90, output_folder='static'):
    """Generates and saves the Heatmap and Network Graph visualizations."""
    timestamp = int(time.time())
    links = {}
    df_vis = df.copy()
    
    if 'Compared_Genome' in df_vis.columns: 
        df_vis.rename(columns={'Compared_Genome': 'Genome2', 'Reference': 'Genome1'}, inplace=True)
        
    all_genomes = pd.concat([df_vis['Genome1'], df_vis['Genome2']]).unique()
    pivot_col_name = "Forward_Similarity(%)" if "Forward_Similarity(%)" in df_vis.columns else "Reverse_Similarity(%)"
    
    # Heatmap Generation
    heatmap_data = df_vis.pivot(index='Genome1', columns='Genome2', values=pivot_col_name).reindex(index=all_genomes, columns=all_genomes)
    heatmap_data_rev = df_vis.pivot(index='Genome2', columns='Genome1', values=pivot_col_name).reindex(index=all_genomes, columns=all_genomes)
    heatmap_data = heatmap_data.fillna(heatmap_data_rev)
    for genome in all_genomes: heatmap_data.loc[genome, genome] = 100.0
    
    if not heatmap_data.empty:
        plt.figure(figsize=(12, 10))
        # Use only genome IDs for ticks when there are many
        if len(all_genomes) > 50:
            sns.heatmap(heatmap_data, annot=False, cmap="viridis", vmin=0, vmax=100, xticklabels=False, yticklabels=False)
        else:
            sns.heatmap(heatmap_data, annot=False, cmap="viridis", vmin=0, vmax=100)
            
        plt.title("Heatmap of Forward Similarity (%)")
        heatmap_filename = f'heatmap_{timestamp}.png'
        plt.savefig(os.path.join(output_folder, heatmap_filename), dpi=300, bbox_inches='tight')
        plt.close()
        links['heatmap'] = os.path.join(output_folder, heatmap_filename)
    
    # Network Graph Generation
    G = nx.Graph()
    df_vis['Max_Jaccard'] = df_vis[['Jaccard_Forward(%)', 'Jaccard_Reverse(%)']].max(axis=1)
    
    for _, row in df_vis.iterrows():
        if row['Max_Jaccard'] >= threshold:
            g1, g2 = row['Genome1'], row['Genome2']
            G.add_edge(g1, g2, weight=row['Max_Jaccard'])
            
    if G.number_of_edges() > 0:
        pos = nx.spring_layout(G, seed=42)
        plt.figure(figsize=(16, 16))
        
        # Draw nodes and edges (simplified visuals for large datasets)
        nx.draw_networkx_nodes(G, pos, node_size=50, node_color="#d62976", alpha=0.8)
        nx.draw_networkx_edges(G, pos, edge_color="gray", alpha=0.3)
        
        # Only draw labels if there are few nodes
        if G.number_of_nodes() < 50:
            nx.draw_networkx_labels(G, pos, font_size=8, font_color="black")
        
        plt.title(f"Genome Similarity Network (Max Jaccard > {threshold}%)")
        network_filename = f'network_{timestamp}.png'
        plt.savefig(os.path.join(output_folder, network_filename), dpi=300, bbox_inches='tight')
        plt.close()
        links['network'] = os.path.join(output_folder, network_filename)
        
    return links
    
# --- FLASK WEB ROUTE ---
@app.route('/analyze', methods=['POST'])
def analyze():
    # ... (Flask web route logic remains the same for small uploads) ...
    # This part handles the interactive web uploads, which is limited.
    try:
        # Simplified for brevity, assume file handling logic is correct
        kmer_size = int(request.form.get('kmer_size', 21))
        threshold = int(request.form.get('threshold', 90))
        
        # NOTE: Full file handling logic must be present here (omitted for brevity)
        uploaded_files = request.files.getlist('all_files')
        saved_paths = []
        for f in uploaded_files:
            path = os.path.join(UPLOAD_FOLDER, f.filename); f.save(path); saved_paths.append(path)
        
        df = compare_all_genomes(saved_paths, k=kmer_size)
        
        if df.empty:
            return jsonify({'status': 'success', 'message': 'Analysis complete, but no data to display.'}), 200
        
        # NOTE: Links are relative paths for web use
        plot_links_abs = generate_and_save_plots(df, threshold=threshold)
        plot_links_web = {k: '/' + v for k, v in plot_links_abs.items()} # Convert absolute paths to web paths
        table_data = df.to_dict(orient='records')
        
        return jsonify({'status': 'success', 'message': 'Analysis complete.', 'table_data': table_data, 'download_links': plot_links_web})
    except Exception as e:
        print(f"An error occurred during web analysis: {e}")
        return jsonify({'status': 'error', 'message': str(e)}), 500

# ----------------------------------------------------------------------
# NEW PERMANENT SECTION: Command Line Interface (CLI) for Big Data
# ----------------------------------------------------------------------

def run_big_analysis_cli(args):
    """
    Function executed when the user runs the script via the command line for big data.
    """
    data_dir = args.directory
    kmer_size = args.kmer_size
    threshold = args.threshold
    
    print("--- SeqNexus Big Data CLI Analysis Started ---")
    print(f"Directory: {data_dir}, K-mer: {kmer_size}, Threshold: {threshold}%")
    
    if not os.path.isdir(data_dir):
        print(f"ERROR: Data directory not found at: {data_dir}")
        return

    # Scan the directory for all FASTA files
    file_paths = glob.glob(os.path.join(data_dir, '*.fasta')) + glob.glob(os.path.join(data_dir, '*.fna'))
    
    if not file_paths:
        print(f"ERROR: No FASTA files found in {data_dir}. Check file extensions (.fasta or .fna).")
        return

    print(f"Found {len(file_paths)} genomes. Running All-vs-All comparison...")
    
    # Run the core comparison function
    df_results = compare_all_genomes(file_paths, k=kmer_size)

    if df_results.empty:
        print("Analysis complete, but no results generated.")
        return

    # Save data table to a file for the user
    output_tsv = os.path.join(STATIC_FOLDER, f'results_k{kmer_size}_t{threshold}_{int(time.time())}.tsv')
    df_results.to_csv(output_tsv, sep='\t', index=False)
    print(f"Results table saved to: {output_tsv}")


    print("Generating high-resolution plots...")
    plot_links = generate_and_save_plots(df_results, threshold=threshold)
    
    print("\n--- CLI Analysis Complete ---")
    print(f"Plots saved to the '{STATIC_FOLDER}' directory:")
    print(f"  Heatmap: {plot_links.get('heatmap', 'N/A')}")
    print(f"  Network: {plot_links.get('network', 'N/A')}")
    print("---------------------------------")


# --- MAIN EXECUTION BLOCK (Controls CLI vs. Web) ---

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="SeqNexus: Alignment-Free Genomic Similarity Tool. Runs a web demo by default, or big analysis via CLI.")
    
    # Define a group for the command-line analysis mode
    cli_group = parser.add_argument_group('Big Data CLI Mode', 'Use these arguments to run analysis directly on a large file directory.')
    cli_group.add_argument('--directory', '-d', type=str, help='[CLI MODE ONLY] Path to the directory containing FASTA files for analysis.')
    cli_group.add_argument('--kmer-size', '-k', type=int, default=21, help='[CLI MODE ONLY] K-mer size for comparison (default: 21).')
    cli_group.add_argument('--threshold', '-t', type=int, default=90, help='[CLI MODE ONLY] Jaccard similarity threshold for network graph creation (default: 90).')
    
    args, unknown = parser.parse_known_args()
    
    # Check if the user supplied the mandatory directory argument for CLI mode
    if args.directory:
        run_big_analysis_cli(args)
    else:
        # Default behavior: Start the Flask web server
        print("Starting Flask web server for interactive demo...")
        print("To run large data analysis, use: python3 app.py --directory /path/to/genomes")
        app.run(debug=True, port=5000)
