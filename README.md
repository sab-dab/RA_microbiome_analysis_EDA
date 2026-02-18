# RA vs HC Gut Microbiome (16S V3–V4) — Exploratory Data Analysis (EDA)

## Project Overview
This repository contains the **exploratory data analysis (EDA) code** for a rheumatoid arthritis (RA) gut microbiome study using a large-scale **16S rRNA V3–V4 amplicon sequencing** dataset.

The goal of this EDA is to evaluate whether the dataset is suitable for downstream microbiome analysis and machine learning, and to check whether **RA vs Healthy Control (HC)** groups show meaningful differences in microbial diversity and community composition.

At the current stage, this repository includes **EDA scripts only** (no modeling code yet).

---

## EDA Objectives
The EDA workflow focuses on three main objectives:

### 1) Dataset Suitability & Integrity
- Verify ASV overlap between the ASV abundance table and taxonomy table  
- Confirm sample and ASV counts

### 2) Data Quality & Sample-Level QC
- Check for missing values
- Compute sequencing depth per sample
- Identify low-depth samples (QC screening)

### 3) Biological Signal & Group Separability
- Alpha diversity:
  - Richness (ASV count)
  - Shannon diversity
- Genus-level summaries:
  - Relative abundance
  - Mean HC vs RA comparisons
  - Top 20 genera visualization
- Beta diversity:
  - Bray–Curtis distance
  - PCoA visualization
  - PERMANOVA and dispersion testing
- ASV prevalence distribution across samples

---

## Dataset (Open Access)
This project uses the open-access dataset from the following publication:

**Li, J., Xu, J., Jin, J., et al. (2025).**  
*A Comprehensive Dataset on Microbiome Dynamics in Rheumatoid Arthritis from a Large-Scale Cohort Study.*  
**Scientific Data**, 12, 232.

The dataset is publicly available from the authors’ official repository (Figshare / journal-linked archive).

---

## Important Note on Data Files This GitHub repository **does not redistribute the dataset files** (such as `.rds` tables).  
It contains **only analysis code**.

To run this code, please download the dataset directly from the original authors’ open-access repository and place the required `.rds` files into your local working directory.

Example input files used in this EDA:
- `1.ASV.profile.rds`
- `1.taxonomy.info.rds`

## Requirements
This project is written in **R** and uses the following key packages:

- `dplyr`
- `tibble`
- `stringr`
- `ggplot2`
- `vegan`
- `tidyr`

---

## How to Run
1. Download the dataset from the official open-access source.
2. Place the `.rds` files in the same folder as the script (or update the file paths).
3. Run the EDA script in RStudio.

---

## Outputs Generated
The EDA script generates:
- Sequencing depth QC plots
- Alpha diversity plots (Richness, Shannon)
- Genus-level relative abundance summaries
- Top 20 genera barplot (HC vs RA)
- PCoA plot (Bray–Curtis)
- PERMANOVA results table (CSV)
- ASV prevalence summary + histogram

---

## Disclaimer
This repository is intended for academic and research use.  
All credit for the dataset belongs to the original authors.  
## Dataset

This project uses the open-access dataset published in:

Li, J., Xu, J., Jin, J., et al. (2025). *A Comprehensive Dataset on Microbiome Dynamics in Rheumatoid Arthritis from a Large-Scale Cohort Study.* Scientific Data, 12, 232.

The dataset is available from the authors' official repository (Figshare):
[https://figshare.com/articles/dataset/Data_for_publication_in_Scientific_Data/27603876]



---
