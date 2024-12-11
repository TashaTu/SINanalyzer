# SINanalyzer

A toolkit for single-sample network construction, network feature analysis, and network visualization.

---

## Introduction

This package provides comprehensive tools for constructing single-sample networks (SSNs), analyzing network features, and visualizing networks. It supports three distinct methods for single-sample network construction:

1. **SiNE**: A novel single-sample network inference method developed by our team for constructing targeted gene networks.  
   Reference: [Chen et al., Briefings in Bioinformatics, 2023](https://doi.org/10.1093/bib/bbad032).
2. **SSN**: Based on *Liu et al., Nucleic Acids Research, 2016*.  
   Reference: [Liu et al., Nucleic Acids Research, 2016](https://doi.org/10.1093/nar/gkw855).
3. **LIONESS**: Adapted from *Kuijjer et al., iScience, 2019*.  
   Reference: [Kuijjer et al., iScience, 2019](https://doi.org/10.1016/j.isci.2019.02.032).

This package is designed to explore individual biological characteristics in complex systems and facilitate personalized network-based studies in diseases.

---

## Installation

To install the package, use the following command:

```bash
pip install SINanalyser
```

Alternatively, clone this repository and install locally:

```bash
git clone https://github.com/username/your-package.git
cd your-package
pip install .
```

## Usage
Here’s an example to get started with **SiNE**
```python
from SINanalyser import integrated_network_construction
import pandas as pd

# Define the file path to your gene expression matrix
gem_path = 'tests/KICH_all.txt'  # Update this path with the location of your GEM file
gem = pd.read_csv(gem_path, sep='\t', header=0, index_col=0)

# Example: specify a subset of genes (or set to None to include all)
gene_df = pd.read_csv('tests/gene_set_100.txt', sep='\t', header=0, index_col=0)
gene_set =  list(gene_df.index)
sample_df = pd.read_csv('tests/sample_list_10.txt', sep='\t', header=0, index_col=0)
sample_list = list(sample_df.index)
outdir = 'tests/SiNE_SIN'  # Specify an output directory for results

# Run SiNE with the specified parameters
integrated_network_construction.network_construction(case_data=gem,method="SiNE",gene_set=gene_set,sample_list=sample_list,outdir=outdir,output_format="npz")
```

### Supported Methods
1. **SiNE**: Specify method="SiNE"
2. **SSN**: Specify method="SSN" (requires control data)
3. **LIONESS**: Specify method="LIONESS"

