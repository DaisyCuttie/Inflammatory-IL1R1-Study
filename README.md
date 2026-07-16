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

Versions in the comments are the ones the analysis was run with. `monocle3` is
**not** required — the pseudotime Rmd reads pre-computed values rather than
recomputing trajectories (see [Known gaps](#known-gaps)).

### Python

Runs on **system Python** (`/usr/bin/python3`, 3.8.10) — not a conda env. The
`scanpy` conda env on this machine is Python 3.10.20 and lacks scvelo; do not
activate it. The default `python3` Jupyter kernel already resolves correctly.

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

> **Pin these.** scvelo 0.2.5 is old and tightly coupled to numpy 1.23 /
> numba 0.56. A plain `pip install scvelo` pulls a modern release that breaks
> these notebooks.

### External tools

- `velocyto` CLI (`velocyto run10x`)
- 10x reference: `refdata-gex-GRCh38-2020-A` (genes GTF)

## Usage

Run in this order. 
### Clustering Analysis
## 1. Load_Single.R
## 2. Positive_Only_Analysis.Rmd
## 3. Postive_CD3_Only.Rmd
### RNA Velocity and Pseudotime Analysis (use scvelo results)
## 1. merge_loompy.ipynb
## 2. scvelo_RNA_Velocity_2025.ipynb
## 3. IL1R1_Pseudotime_2025.ipynb
## 4. pseudotime_modulescore_plotting.Rmd
