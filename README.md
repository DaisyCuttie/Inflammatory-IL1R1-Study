# Inflammatory IL1R1 Study

Analysis of IL1R1+ signaling in vivo 10x scRNA seq experiements

## Overview

Scripts of loading 10x scRNA seq data, downstream analysis, marker lookup, RNA velocity analysis and pseudotime anlysis

## Requirements

### R

R 4.3.2. Install from CRAN:

```r
install.packages(c(
  "Seurat",      # 5.4.0.9022
  "dplyr",       # 1.1.4
  "ggplot2",     # 4.0.2
  "harmony",     # 1.2.4
  "clustree",    # 0.5.1
  "readr",       # 2.1.6
  "openxlsx",    # 4.2.8.1
  "tidyr",       # 1.3.1
  "scales",      # 1.4.0
  "Matrix",      # 1.6.5
  "knitr", "rmarkdown"
))
```

### Python

Runs on Python 3.8.10.

```bash
pip install \
  scvelo==0.2.5 \
  velocyto==0.17.17 \
  loompy==3.0.7 \
  scanpy==1.9.1 \
  anndata==0.8.0 \
  matplotlib==3.5.2 \
  numpy==1.23.5 \
  pandas==1.4.2 \
  numba==0.56.4 \
  llvmlite==0.39.1 \
  h5py==3.7.0
```

### External tools

- `velocyto` CLI (`velocyto run10x`)
- 10x reference: `refdata-gex-GRCh38-2020-A` (genes GTF)

## Usage

Run in this order. 
### Clustering Analysis
```bash
Load_Single.R	# load in 10x cellranger outputs and do some standard scRNA analysis
```
```bash
2. Positive_Only_Analysis.Rmd	# some exploratory studies and re-analyze positive data only
```
```bash
3. Postive_CD3_Only.Rmd	# some other exploratory studies and made decision with using positive CD3 only
```
### RNA Velocity and Pseudotime Analysis (use scvelo results)
```bash
1. merge_loompy.ipynb	# process loom file outputing from the scvelo and merge all samples
```
```bash
2. scvelo_RNA_Velocity_2025.ipynb	# subsetting cells that are only in the seurat object
```
```bash
3. IL1R1_Pseudotime_2025.ipynb	# pseudotime analysis based on the RNA velocity
```
```bash
4. pseudotime_modulescore_plotting.Rmd	# module score plotting along the pseudotime trend
```
