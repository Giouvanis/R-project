RNA-Seq Analysis: Weathering Paper
This repository contains the transcriptomics pipeline and principal component analysis (PCA) for my research on biological weathering. The study explores the differential gene expression and variance between symbiotic and non-symbiotic states across different mineral conditions.

📊 Analysis Overview
The analysis is performed in R using the DESeq2 framework. Key components include:

Differential Expression: Utilizing a 2-factorial design (Mineral + Symbiotic + Mineral:Symbiotic) to identify significant genetic responses.

3D Visualization: Implementation of interactive 3D PCA plots using plotly and scatterplot3d to visualize sample clustering.

Data Normalization: Variance stabilizing transformations (VST) and regularized log transformation (rlog) for downstream visualization.

🛠 Tech Stack
Language: R

Bioinformatics: Bioconductor, DESeq2, edgeR

Visualization: ggplot2, plotly, scatterplot3d, RColorBrewer

📂 Project Structure
RNA_seq_weathering_paper_PCA_analysis.R: The main analysis script containing the DESeq2 pipeline and plotting functions.

raw_counts.csv: The transcriptomic count matrix.

sample_list.csv: Metadata containing the experimental design (Symbiotic vs. Non-Symbiotic; Mineral Plus vs. Minus).

PCA.csv: Exported principal component data for reproducibility.

📈 Key Results
The analysis identifies how the presence of specific minerals and symbiotic relationships drive global transcriptomic changes.

(Pro-tip: After you upload your PCA plot image to GitHub, you can drag and drop the image right here in the README to show the results visually!)
