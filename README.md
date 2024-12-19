# SINanalyzer

A toolkit for single-sample network construction, network feature analysis, and network visualization.

---

## Introduction

This package provides comprehensive tools for constructing single-sample networks, analyzing network features, and visualizing networks. It supports three distinct methods for single-sample network construction:

1. **SiNE**: A novel single-sample network inference method developed by our team for constructing targeted gene networks.  
   Reference: [Chen et al., Briefings in Bioinformatics, 2023](https://doi.org/10.1093/bib/bbad032).
2. **SSN**: Based on *Liu et al., Nucleic Acids Research, 2016*.  
   Reference: [Liu et al., Nucleic Acids Research, 2016](https://doi.org/10.1093/nar/gkw855).
3. **LIONESS**: Adapted from *Kuijjer et al., iScience, 2019*.  
   Reference: [Kuijjer et al., iScience, 2019](https://doi.org/10.1016/j.isci.2019.02.032).

This package is designed to explore individual biological characteristics in complex systems and facilitate personalized network-based studies in diseases.

---

## Installation

1. Use the following command to clone this repository and install locally

```bash
git clone https://github.com/TashaTu/SINanalyzer.git
cd SINanalyzer
```
2. Create a virtual environment (optional but recommended)
```bash
python3 -m venv env
source env/bin/activate  # For MacOS/Linux
# OR
.\env\Scripts\activate   # For Windows
```
3. Install the package
```
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

# Example: Specify gene set and sample list (optional)
gene_df = pd.read_csv('tests/gene_set_100.txt', sep='\t', header=0, index_col=0) # The example file format with a header, so you should set `header=0`.
gene_set =  list(gene_df.index)
sample_df = pd.read_csv('tests/sample_list_10.txt', sep='\t', header=0, index_col=0) # The example file format with a header, so you should set `header=0`.
sample_list = list(sample_df.index)
outdir = 'tests/SiNE_SIN'  # Specify an output directory for results

# Run SiNE with the specified parameters
integrated_network_construction.network_construction(case_data=gem,method="SiNE",gene_set=gene_set,sample_list=sample_list,outdir=outdir,output_format="npz")
```
## Parameters

### Core Parameters

- **`case_data` (pd.DataFrame)**:  
  A DataFrame containing the case gene expression data. The data should be formatted with genes as rows and samples as columns.  
  **Required**.

- **`control_data` (pd.DataFrame or None)**:  
  A DataFrame containing the control gene expression data. Only required for the SSN method. If not used, the default is `None`.

- **`method` (str)**:  
  The network construction method to use. Options are:
  - `"SiNE"`
  - `"LIONESS"`
  - `"SSN"`

- **`amp` (float)**:  
  Amplitude parameter for SiNE network construction. Default is `0.1`.

- **`pvalue` (float)**:  
  P-value threshold for network construction. Default is `0.05`.

### Additional Parameters

- **`gene_set` (list or None)**:  
  A list of genes to focus on. If `None`, all genes in the dataset are used.

- **`sample_list` (list or None)**:  
  A list of sample names to include in the output. If `None`, all samples are output.

- **`normalize` (bool)**:  
  If `True`, the data will be log2 normalized. Default is `False`.

- **`outdir` (str or None)**:  
  Directory to save the output files. If `None`, the current directory is used.

- **`output_format` (str)**:  
  Format for output files. Options are:
  - `"npz"`
  - `"edge_list"`
  - `"edge_list_score"`
  - `"edge_list_zscore"`

### Notes

- This function dynamically adjusts parameters and execution based on the specified method.  
- Make sure to provide all required inputs, especially when using methods like **SSN**, which require both `case_data` and `control_data`.

### More Details
For more detailed usage and examples, please refer to the `Package_usage_guide.ipynb` file.

