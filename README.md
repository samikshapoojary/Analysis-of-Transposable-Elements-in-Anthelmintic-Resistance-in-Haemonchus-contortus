# Exploring the Role of Transposable Elements in Anthelmintic Resistance in *Haemonchus contortus*

**MSc Bioinformatics Project – University of Glasgow**  
**Author:** Samiksha Poojary  
**Supervisors:** Dr Roz Laing & Prof James Cotton  

---

## 📘 Overview

This repository contains all scripts developed for the MSc Bioinformatics project *“Exploring the Role of Transposable Elements in Anthelmintic Resistance in Haemonchus contortus.”*  
All code is provided for reproducibility and as reference material.

The project investigates transposable element (TE) dynamics and their potential involvement in anthelmintic resistance, integrating results from **EarlGrey** and **PoPoolationTE2**, with supplementary RNA-seq workflows included in the appendix.

---

## 📂 Repository Contents

### **A. Master Scripts**

#### `master_bash_script.sh`
Contains all Bash workflows used for:
- FASTQ header/format fixing  
- **EarlGrey** TE library generation  
- Filtering TE annotations  
- Preparing merged genome + TE FASTA and TE hierarchy file  
- **PoPoolationTE2** workflow (sample-level and joint analysis)  
- Overlap analysis with genes and UTRs  

#### `master_r_script.R`
Contains all R workflows used for:
- Parsing and cleaning TE annotation data  
- Summarising TE content across genomes  
- Plotting genome composition, TE landscapes, and distributions  
- Statistical analyses (logistic regression, PCA, UpSet plots, etc.)  

---

### **A.2. RNA-seq Scripts (Supplementary)**

#### `supp_rnaseq_bash.sh`
Contains RNA-seq analysis steps for:
- **STAR** genome index building  
- **RNA-seq read alignment**  
- **HTSeq-count** for read quantification  

> **Note:** RNA-seq steps were not part of the main thesis analysis but are included here as supplementary material for completeness.

---

## ⚙️ Requirements

### **B.1. Software**
- Bash environment (SLURM scheduler used for job submission)
- R with the following libraries:
  - `tidyverse` (dplyr, ggplot2, tibble, tidyr, readr, forcats)
  - `GenomicRanges`, `IRanges`
  - `ComplexUpset`, `patchwork`, `ggpubr`, `scales`, `RColorBrewer`

### **B.2. Bioinformatics Tools**
- [EarlGrey](https://github.com/TobyBaril/EarlGrey/blob/main/README.md) 
- [PoPoolationTE2](https://sourceforge.net/p/popoolation-te2/wiki/Home/)  
- [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)  
- [BEDTools](https://bedtools.readthedocs.io/en/latest/)  
- [SAMtools](https://www.htslib.org/doc/samtools.html)  
- [STAR](https://github.com/alexdobin/STAR)  
- [HTSeq-count](https://d.docs.live.net/A7A1307BF1E6342B/Documents/htseq-count%20manual)  
---

## 🧬 Input Data

The project used the following genome assemblies and sequence datasets:

| Organism | Assembly ID | Source |
|-----------|--------------|--------|
| *Haemonchus contortus* | PRJEB506 (WBPS18) | WormBase ParaSite |
| *Teladorsagia circumcincta* | WSI3.0 | WormBase ParaSite |
| *Caenorhabditis elegans* | PRJNA13758 (WBPS19) | WormBase ParaSite |

**Sequencing data:**
- 9 pooled genomic samples (Pool-seq)  
- RNA-seq FASTQ files (not analysed in final thesis)
- Pooled whole-genome sequencing (Pool-seq) data and RNA-seq datasets for control, ivermectin-resistant, and moxidectin-resistant Haemonchus contortus populations were obtained from the Centre for Genomic Research (CGR), University of Liverpool:
[Pool-seq](https://cgr.liv.ac.uk/illum/LIMS28593_fbc6c22b786c7600/) 
[RNA-seq](https://cgr.liv.ac.uk/illum/LIMS30731_cb27564db3511460/)

---

## ▶️ How to Run

### **Genome & TE Analysis (Bash)**
Run `master_bash_script.sh`, executing each section corresponding to a specific analysis step (e.g., genome preparation, PoPoolationTE2 runs, feature overlap analysis).

### **R Analysis & Plotting**
Run `master_r_script.R` section by section in R or RStudio for data parsing, visualization, and statistical analysis.

> ⚠️ File paths within scripts are project-specific (e.g. `/Data3/samiksha/...`). Adjust paths according to your own environment before running.

---

## 🧾 Notes

- The master scripts are concatenated versions of all smaller scripts used in the project for convenience and assessment.  
- It is recommended to execute each section separately rather than running the entire master file at once.  
- Data files (e.g. FASTQ, BAM, and intermediate results) are **not included** due to size constraints.


GitHub: [https://github.com/YOURUSERNAME](https://github.com/YOURUSERNAME)

