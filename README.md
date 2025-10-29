# **🧬 SeqNexus: Alignment-Free Genomic Data Analysis Platform**

[![GitHub language count](https://img.shields.io/github/languages/count/divya0002/seqnexus-app)](https://github.com/divya0002/seqnexus-app)
[![GitHub top language](https://img.shields.io/github/languages/top/divya0002/seqnexus-app)](https://github.com/divya0002/seqnexus-app)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

---

## **💡 Project Overview: Bridging Genomics and Data Science**

**SeqNexus** is an interactive web tool built on Python and Flask for **high-throughput, alignment-free comparison of genomic sequences**. It utilises k-mer-based metrics to rapidly derive genetic similarity, providing a crucial exploratory data analysis (EDA) tool for fields like microbial genomics and epidemiology (e.g., analysing strains like *M. tuberculosis*).

This tool is a demonstration of proficiency in:

* **Bioinformatics:** Handling, parsing, and computationally analyzing large FASTA files.
* **Data Science:** Feature engineering (k-mer profiles), data modeling, and visualization for complex biological relationships.
* **Web Development:** Implementing a functional user interface and backend using a Flask RESTful API.

## **📈 Data Science & AI/ML Vision: Universal MTB Genome Predictor**

The ultimate goal of SeqNexus is to evolve into a **Universal Tool for *M. tuberculosis* Genome Prediction**, leveraging the extracted k-mer/SSR features alongside industry-standard **Single Nucleotide Polymorphism (SNP)** features derived from public catalogs (like WHO/CRYPTIC).

This multi-model architecture focuses on three high-impact predictions:

### 1. Action 1: Universal Lineage Predictor

* **Goal:** Classify the uploaded genome into one of the **9 universally accepted MTB Lineages (L1–L9)**.
* **Features:** Primarily uses **Lineage-Defining SNPs** (LSP-defined mutations) and secondarily uses your **SSR counts** to fine-tune prediction within a lineage.
* **Model:** A multi-class classifier (e.g., **Random Forest**) trained on a large public dataset (CRYPTIC or TB Portals).

### 2. Action 2: Geographic Origin Predictor

* **Goal:** Predict the **continent or country of origin** of the strain.
* **Features:** Exploits known **Lineage-to-Geography Correlations** and uses the predicted **MTB Lineage** (from Action 1) combined with your **SSR Counts** (as SSRs can help define geographically restricted sub-clusters).
* **Model:** A classifier trained on public metadata linking genome accession numbers to the Country or WHO Region of isolation.

### 3. Action 3: Universal Drug Resistance Predictor

* **Goal:** Predict resistance for common anti-TB drugs (e.g., Isoniazid, Rifampicin).
* **Features:** Primarily uses the **Binary presence/absence of WHO-defined resistance-conferring mutations** (e.g., in genes like *rpoB*, *katG*, *pncA*). Secondary features include your **SSR features in resistance hotspot genes**.
* **Model:** 13 separate **Binary Classifiers** (one for each drug) trained on the massive CRYPTIC dataset.

---

## **✨ Core Features**

* **Alignment-Free Metrics:** Efficient calculation of **Jaccard Index** and **Raw Similarity** based on k-mer co-occurrence.
* **Batch Analysis Modes:** Supports two primary research workflows:
    * **All-vs-All:** Comparing every genome against every other genome in the uploaded cohort.
    * **Reference-vs-Query:** Comparing a single known reference sequence against multiple unknown query sequences.
* **User-Defined Parameters:** Allows dynamic adjustment of **k-mer size** and the **similarity threshold** for network visualization.

## **🖼️ Demo Visualizations**

The output of SeqNexus is designed to facilitate rapid analytical insight through two key graphics, downloadable in high-resolution. *(NOTE: You should replace the placeholder lines below with actual image links after you generate and upload your demo files to the repository.)*

**Figure 1: Similarity Heatmap**
*A grid visualization of the pairwise genetic similarity (e.g., Forward Similarity %) between every genome in the analysis set.*
![Similarity Heatmap Example Image Placeholder](assets/heatmap_example.png)

**Figure 2: Genome Relationship Network**
*A network graph showing relationships between genomes that meet a user-defined similarity threshold, visualizing potential clusters or epidemiological links.*
![Genome Network Graph Example Image Placeholder](assets/network_graph_example.png)

## **⚙️ Technical Stack**

| Component | Technology | Role in Project | 
| :--- | :--- | :--- | 
| **Backend Framework** | `Python`, `Flask` | Provides the lightweight web server and RESTful API endpoint (`/analyze`). | 
| **Sequence Handling** | `Biopython` (SeqIO, Seq) | Robust tools for reading FASTA files, sequence manipulation, and reverse complement calculation. | 
| **Data Processing** | `Pandas`, `Collections` | Efficient tabulation of results and k-mer counting/frequency analysis. | 
| **Visualization** | `Matplotlib`, `Seaborn`, `NetworkX` | Libraries used to generate static, publication-ready images (Heatmaps and Graphs). | 

## **🚀 Local Installation and Usage**

To run SeqNexus locally and avoid memory/scaling issues associated with free cloud tiers (especially for large genomes), follow these steps:

1.  **Clone the Repository:**

    ```bash
    git clone [https://github.com/divya0002/seqnexus-app.git](https://github.com/divya0002/seqnexus-app.git)
    cd seqnexus-app
    ```

2.  **Install Dependencies:**

    ```bash
    # Ensure you have Python 3.9+ installed
    pip install -r requirements.txt
    ```

3.  **Run the Web Server:**

    ```bash
    python app.py
    ```

    *The application will now be running on your local machine. Open your web browser and navigate to `http://127.0.0.1:5000/` to use the interactive interface.*

## **⚠️ Project Status & Scalability Note**

SeqNexus is maintained as a **Proof-of-Concept** with all its code available in the repository. The currently deployed web application is limited to **small-to-medium datasets** for demonstration purposes due to hosting memory constraints. The core Python architecture is designed for scalability and can be adapted for parallel processing and deployment on High-Performance Computing (HPC) environments to handle full-sized genomes.
