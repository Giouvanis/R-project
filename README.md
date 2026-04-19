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


<img width="972" height="715" alt="Picture_1 (PC1 and PC2 variance)" src="https://github.com/user-attachments/assets/810dde66-1bd0-4442-a56d-dddac800600d" />
<img width="1065" height="643" alt="Picture_3 group_2" src="https://github.com/user-attachments/assets/9ace30c0-c3e5-4f55-a424-1ce5c086e294" />
<img width="984" height="624" alt="Picture_2 group" src="https://github.com/user-attachments/assets/b8a0e591-c815-473e-b019-cc644b0555bf" />


