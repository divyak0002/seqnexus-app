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

app = Flask(__name__)
CORS(app)

UPLOAD_FOLDER = 'uploads'
STATIC_FOLDER = 'static'
# ## THIS IS THE FIX: Added exist_ok=True to prevent race conditions ##
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(STATIC_FOLDER, exist_ok=True)


@app.route('/')
def home():
    return render_template('i.html')

def load_sequence(fasta_path):
    with open(fasta_path, "r") as handle:
        record = next(SeqIO.parse(handle, "fasta"))
        return str(record.seq).upper()

def generate_kmers(sequence, k):
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

def compute_similarity(seq1, seq2, k=21):
    kmers1, kmers2 = generate_kmers(seq1, k), generate_kmers(seq2, k)
    if not kmers1 or not kmers2: return 0, 0
    counts1, counts2 = Counter(kmers1), Counter(kmers2)
    common_kmers = set(counts1) & set(counts2)
    shared = sum(min(counts1[k_val], counts2[k_val]) for k_val in common_kmers)
    forward_sim = (shared / len(kmers1)) * 100 if len(kmers1) > 0 else 0
    union_kmers = set(counts1) | set(counts2)
    jaccard_sim = (len(common_kmers) / len(union_kmers)) * 100 if len(union_kmers) > 0 else 0
    return forward_sim, jaccard_sim

def compare_reference_vs_all(reference_path, query_paths, k=21):
    ref_seq, ref_id = load_sequence(reference_path), os.path.basename(reference_path)
    results = []
    for query_path in query_paths:
        query_id, query_seq = os.path.basename(query_path), load_sequence(query_path)
        query_seq_rc = str(Seq(query_seq).reverse_complement())
        forward_sim, jaccard_sim = compute_similarity(ref_seq, query_seq, k)
        reverse_sim, jaccard_sim_rc = compute_similarity(ref_seq, query_seq_rc, k)
        results.append([ref_id, query_id, round(forward_sim, 2), round(reverse_sim, 2), round(jaccard_sim, 2), round(jaccard_sim_rc, 2)])
    df = pd.DataFrame(results, columns=["Reference", "Compared_Genome", "Forward_Similarity(%)", "Reverse_Similarity(%)", "Jaccard_Forward(%)", "Jaccard_Reverse(%)"])
    df = df[["Reference", "Compared_Genome", "Jaccard_Forward(%)", "Jaccard_Reverse(%)", "Reverse_Similarity(%)", "Forward_Similarity(%)"]]
    return df

def compare_all_genomes(file_paths, k=21):
    sequences = {os.path.basename(p): load_sequence(p) for p in file_paths}
    filenames = list(sequences.keys())
    results = []
    for i in range(len(filenames)):
        for j in range(i + 1, len(filenames)):
            file1, file2 = filenames[i], filenames[j]
            seq1, seq2 = sequences[file1], sequences[file2]
            seq2_rc = str(Seq(seq2).reverse_complement())
            forward_sim, jaccard_sim = compute_similarity(seq1, seq2, k)
            reverse_sim, jaccard_sim_rc = compute_similarity(seq1, seq2_rc, k)
            results.append([file1, file2, round(forward_sim, 2), round(reverse_sim, 2), round(jaccard_sim, 2), round(jaccard_sim_rc, 2)])
    df = pd.DataFrame(results, columns=["Genome1", "Genome2", "Forward_Similarity(%)", "Reverse_Similarity(%)", "Jaccard_Forward(%)", "Jaccard_Reverse(%)"])
    df = df[["Genome1", "Genome2", "Jaccard_Forward(%)", "Jaccard_Reverse(%)", "Reverse_Similarity(%)", "Forward_Similarity(%)"]]
    return df

def generate_and_save_plots(df, threshold=80, output_folder='static'):
    timestamp = int(time.time())
    links = {}
    df_vis = df.copy()
    if 'Compared_Genome' in df_vis.columns: df_vis.rename(columns={'Compared_Genome': 'Genome2', 'Reference': 'Genome1'}, inplace=True)
    all_genomes = pd.concat([df_vis['Genome1'], df_vis['Genome2']]).unique()
    pivot_col_name = "Forward_Similarity(%)" if "Forward_Similarity(%)" in df_vis.columns else "Reverse_Similarity(%)"
    heatmap_data = df_vis.pivot(index='Genome1', columns='Genome2', values=pivot_col_name).reindex(index=all_genomes, columns=all_genomes)
    heatmap_data_rev = df_vis.pivot(index='Genome2', columns='Genome1', values=pivot_col_name).reindex(index=all_genomes, columns=all_genomes)
    heatmap_data = heatmap_data.fillna(heatmap_data_rev)
    for genome in all_genomes: heatmap_data.loc[genome, genome] = 100.0
    if not heatmap_data.empty:
        plt.figure(figsize=(10, 8))
        sns.heatmap(heatmap_data, annot=False, cmap="viridis", vmin=0, vmax=100)
        plt.title("Heatmap of Forward Similarity (%)")
        heatmap_filename = f'heatmap_{timestamp}.png'
        plt.savefig(os.path.join(output_folder, heatmap_filename), dpi=150, bbox_inches='tight')
        plt.close()
        links['heatmap'] = f'/static/{heatmap_filename}'
    G = nx.Graph()
    for _, row in df.iterrows():
        if row[pivot_col_name] >= threshold:
            g1, g2 = row.get("Genome1", row.get("Reference")), row.get("Genome2", row.get("Compared_Genome"))
            G.add_edge(g1, g2, weight=row[pivot_col_name])
    if G.number_of_edges() > 0:
        pos = nx.spring_layout(G, seed=42)
        plt.figure(figsize=(12, 12))
        nx.draw(G, pos, with_labels=True, node_size=1000, node_color="#d62976", edge_color="gray", font_size=10, font_color="white")
        labels = nx.get_edge_attributes(G, 'weight')
        nx.draw_networkx_edge_labels(G, pos, edge_labels={k: f"{v:.1f}%" for k, v in labels.items()})
        plt.title(f"Genome Similarity Network (threshold > {threshold}%)")
        network_filename = f'network_{timestamp}.png'
        plt.savefig(os.path.join(output_folder, network_filename), dpi=150, bbox_inches='tight')
        plt.close()
        links['network'] = f'/static/{network_filename}'
    return links
    
@app.route('/analyze', methods=['POST'])
def analyze():
    try:
        mode = request.form.get('mode', 'All vs All')
        kmer_size = int(request.form.get('kmer_size', 21))
        threshold = int(request.form.get('threshold', 90))
        uploaded_files = request.files.getlist('all_files') if mode == 'All vs All' else request.files.getlist('query_files')
        reference_file = request.files.get('reference_file') if mode == 'Pairwise' else None
        
        # Clear the uploads folder before saving new files
        for f in os.listdir(UPLOAD_FOLDER):
            os.remove(os.path.join(UPLOAD_FOLDER, f))

        saved_paths = []
        for f in uploaded_files:
            path = os.path.join(UPLOAD_FOLDER, f.filename); f.save(path); saved_paths.append(path)
        reference_path = None
        if reference_file:
            reference_path = os.path.join(UPLOAD_FOLDER, reference_file.filename); reference_file.save(reference_path)
        
        df = compare_all_genomes(saved_paths, k=kmer_size) if mode == 'All vs All' else compare_reference_vs_all(reference_path, saved_paths, k=kmer_size)
        
        if df.empty:
            return jsonify({'status': 'success', 'message': 'Analysis complete, but no data to display.'}), 200
        
        plot_links = generate_and_save_plots(df, threshold=threshold)
        table_data = df.to_dict(orient='records')
        return jsonify({'status': 'success', 'message': 'Analysis complete.', 'table_data': table_data, 'download_links': plot_links})
    except Exception as e:
        # For debugging, let's print the exception to the log
        print(f"An error occurred during analysis: {e}")
        return jsonify({'status': 'error', 'message': str(e)}), 500

# This part is removed for deployment as Gunicorn starts the app
# if __name__ == '__main__':
#     app.run(debug=True, port=5001)
