# 🧬 **SeqNexus: Alignment-Free Genomic Data Analysis Platform**

---

## 💡 **Project Overview: Bridging Genomics and Data Science**

SeqNexus is an interactive web tool built on **Python** and **Flask** for **high-throughput, alignment-free comparison** of genomic sequences. It utilises **k-mer-based metrics** to rapidly derive genetic similarity, providing a crucial **Exploratory Data Analysis (EDA)** tool for fields like microbial genomics and epidemiology (e.g., analysing strains like *M. tuberculosis*).

This tool is a demonstration of proficiency in:

* **Bioinformatics:** Handling, parsing, and computationally analyzing large FASTA files.
* **Data Science:** Feature engineering (k-mer profiles), data modeling, and visualization for complex biological relationships.
* **Web Development:** Implementing a functional user interface and backend using a Flask RESTful API.

---

## 📈 **Data Science & AI/ML Vision: Universal MTB Genome Predictor**

The ultimate goal of SeqNexus is to evolve into a **Universal Tool for *M. tuberculosis* Genome Prediction**, leveraging the extracted k-mer/SSR features alongside industry-standard Single Nucleotide Polymorphism (SNP) features derived from public catalogs (like WHO/CRYPTIC).

This multi-model architecture focuses on three high-impact predictions:

### 1. **Action 1: Universal Lineage Predictor**
* **Goal:** Classify the uploaded genome into one of the 9 universally accepted MTB Lineages (L1–L9).
* **Features:** Primarily uses Lineage-Defining SNPs (LSP-defined mutations) and secondarily uses your SSR counts to fine-tune prediction within a lineage.
* **Model:** A multi-class classifier (e.g., Random Forest) trained on a large public dataset (CRYPTIC or TB Portals).

### 2. **Action 2: Geographic Origin Predictor**
* **Goal:** Predict the continent or country of origin of the strain.
* **Features:** Exploits known Lineage-to-Geography Correlations and uses the predicted MTB Lineage (from Action 1) combined with your SSR Counts (as SSRs can help define geographically restricted sub-clusters).
* **Model:** A classifier trained on public metadata linking genome accession numbers to the Country or WHO Region of isolation.

### 3. **Action 3: Universal Drug Resistance Predictor**
* **Goal:** Predict resistance for common anti-TB drugs (e.g., Isoniazid, Rifampicin).
* **Features:** Primarily uses the Binary presence/absence of WHO-defined resistance-conferring mutations (e.g., in genes like *rpoB*, *katG*, *pncA*). Secondary features include your SSR features in resistance hotspot genes.
* **Model:** 13 separate Binary Classifiers (one for each drug) trained on the massive CRYPTIC dataset.

---

## ✨ **Core Features**

* **Alignment-Free Metrics:** Efficient calculation of **Jaccard Index** and **Raw Similarity** based on k-mer co-occurrence.
* **Batch Analysis Modes:** Supports two primary research workflows:
    * **All-vs-All:** Comparing every genome against every other genome in the uploaded cohort.
    * **Reference-vs-Query:** Comparing a single known reference sequence against multiple unknown query sequences.
* **User-Defined Parameters:** Allows dynamic adjustment of **k-mer size** and the **similarity threshold** for network visualization.

---

## 🖼️ **Strategically Chosen Visualizations**

We provide two images to showcase different strengths. These images were generated using optimized genome snippets for speed and visual clarity.

* **Asset 1:** Scale and Robustness (100 Genomes)

This Heatmap proves the tool can process and correctly map 100 and even more simultaneous comparisons, demonstrating computational scale.

| Image Asset | Source Data | Proof |
| :--- | :--- | :--- |
| **Similarity Heatmap** | 100 Genomes (Snippets) | Scale and Computational Robustness |

* **Asset 2:** Clustering and Clarity (5 Genomes)

This Network Graph uses fewer nodes to clearly illustrate the maximum Jaccard similarity threshold, revealing distinct clusters and making visual analysis possible.

| Image Asset | Source Data | Proof |
| :--- | :--- | :--- |
| **Network Graph** | 5 Genomes (Snippets) | Visualization Clarity and Clustering Logic |

---

## ⚙️ **Technical Stack**

| Component | Technology | Role in Project |
| :--- | :--- | :--- |
| **Backend Framework** | Python, Flask | Provides the lightweight web server and RESTful API endpoint (`/analyze`). |
| **Sequence Handling** | Biopython (SeqIO, Seq) | Robust tools for reading FASTA files, sequence manipulation, and reverse complement calculation. |
| **Data Processing** | Pandas, Collections | Efficient tabulation of results and k-mer counting/frequency analysis. |
| **Visualization** | Matplotlib, Seaborn, NetworkX | Libraries used to generate static, publication-ready images (Heatmaps and Graphs). |

---

## 🚀 **Local Installation and Usage**

The SeqNexus tool supports two distinct execution modes: an interactive web demo and a dedicated Command-Line Interface (CLI) for large datasets.

### 1. **Setup (Mandatory)**

Clone the repository and install all required libraries:

```bash
# Clone the repository (Use SSH for reliability if configured)
git clone git@github.com:divya0002/seqnexus-app.git
# Navigate and Install Dependencies
cd seqnexus-app
pip install -r requirements.txt
```

### 2. **Mode A: Interactive Web Demo (Small Data)**

Use this mode for testing and visualizing small datasets directly in your browser.

```bash
# Start the Flask web server
python3 app.py
```
* **Public Access:** The live demo can be accessed here: https://seqnexus-app.onrender.com

* **Local Access (Developer only):** Open your browser and navigate to http://127.0.0.1:5000/.

* **Use Case:** Interactive testing, quick comparisons ($\le 10$ files), and UI demonstration.

### 3. **Mode B: Production Scale Analysis (CLI)**

Use this mode for processing large genome cohorts (e.g., 400+ genomes). This mode bypasses the web interface and saves the results directly to the disk, which is more stable for heavy computation.

### Demonstration Commands

```bash
# Command for 100-Genome Heatmap (Scale Proof)
# NOTE: This uses the default threshold of 95% (as defined in app.py)
python3 app.py -d [path/to/100_genome_snippets]

# Command for 5-Genome Network (Clarity Proof)
# NOTE: Use a smaller dedicated folder with only 5 files for clearer visualization.
python3 app.py -d [path/to/5_genome_snippets] -t 95
```

### CLI Arguments

| Argument | Description | Default |
| :--- | :--- | :--- |
| ```bash #--directory / -d ``` | REQUIRED. Path to the folder containing all FASTA files. | N/A |
| ```bash #--mode / -m ``` | Comparison mode: `/all-vs-all` or `/pairwise` (Reference-vs-All). | `/all-vs-all` |
| ```bash #--reference-file / -r ``` | $$Pairwise Mode Only:$$ Name of the specific reference file within the --directory to compare others against. | First file found |
| ```bash #--kmer-size / -k ``` | Length of the k-mers (features) to use for comparison. | 21 |
| ```bash #--threshold / -t ``` | Minimum Max-Jaccard similarity (%) required to draw a link in the network graph. Use 95+ for large datasets to reduce clutter. | 90 (or 95 in code) |

# **⏱️ Execution Time Note**

Large-scale All-vs-All comparison is computationally intensive (complexity $\propto N^2$).
* **Example Benchmark:** Analyzing $\sim$400 genomes can take 30 to 60+ minutes on a standard desktop CPU, as it involves over 95,000 similarity calculations. The script is running normally if it produces no output for extended periods.

# **⚠️ Project Status & Scalability Note**

SeqNexus is maintained as a Proof-of-Concept with all its code available in the repository. The currently deployed web application is limited to small-to-medium datasets for demonstration purposes due to hosting memory constraints. The core Python architecture is designed for scalability and can be adapted for parallel processing and deployment on High-Performance Computing (HPC) environments to handle full-sized genomes.
