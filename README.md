# scEpiDel

An algorithm for the detection of full genomic deletions in single-cell ATAC-seq cancer datasets. It identifies large-scale copy number losses by comparing read distributions between cancer and normal cell populations using a variable-window size that adapts to its local read density. Developed by the Perkins Lab at the Ottawa Hospital Research Institute.

## Overview

Single-cell ATAC-seq (scATAC-seq) enables chromatin accessibility profiling at single-cell resolution.  
However, detecting full genomic deletions remains challenging due to sparse read coverage and uneven signal distribution across cells.

`scEpiDel` addresses this by:

- Aggregating fragment data by cell type (cancer vs normal)
- Removing overlapping reads
- Sliding a window that adapts to local read density across the genome
- Computing Z-scores on read distributions to detect regions significantly depleted in cancer cells
- Merging adjacent significant windows to identify continuous deletion regions

### Description of the algorithm:

    1.  Split fragment file into separate bed files for reads from cancer cells and non-cancer cells
    2.  Remove overlapping reads from both files (if a read began within the genomic interval of the previous read, it was discarded to ensure that only non-overlapping fragments were retained.)
    3.  Scan the genome using a sliding window containing N consecutive cancer reads, advancing the window by K reads at each step. Count the number of non-cancer reads within the window
    4.  Apply a logarithmic transformation to the non-cancer read counts and standardize them across all windows to compute Z-scores
    5.  Identify windows with Z-scores exceeding a defined significance threshold
    6.  Merge consecutive significant windows to form continuous deletion regions
    7.  Output the genomic coordinates and statistical summaries of the predicted deletion intervals

## Installation

```bash
pip install scEpiDel
```

## Usage

### Example Usage:

```python
from scEpiDel import sort_fragments, scEpiDel

chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

# Sort fragments by cell type and remove overlapping reads. Returns paths to filtered fragment files.
cancer_file, normal_file = sort_fragments(fragmentFile = "data/fragments.tsv", # path to fragment file
                                          cellTypeFile = "data/celltypes.csv", # path to cell type file with barcodes and cancer cells labeled
                                          cancerLabel = "Cancer", # label used for cancer cells in cell type file
                                          outputDir = "output/", # directory to save filtered fragment files
                                          chromosomes = chromosomes) # list of chromosomes to process

# Run scEpiDel. Returns dataframe with predicted deletions and statisitcs.
results = scEpiDel(cancer_file = cancer_file, # path to filtered cancer fragment file
                   normal_file = normal_file, # path to filtered normal fragment file
                   chromosomes = chromosomes, # list of chromosomes to analyze
                   window_size = 10, # size of sliding window in # of cancer cell reads
                   z_cutoff = 4.5) # Z-score cutoff for significant deletions

# Save results
results.to_csv("output/deletions.tsv", sep="\t", index=False)
```
### Input Formats

Fragment File: 
```bash
#Standard fragment file format

$ head fragments.tsv 
chr1    10001   10101   TAGTACGGTTAACAGT-1      1
chr1    10006   10204   TGGTGATTCTCAATTC-1      1
chr1    10006   10386   ATGAAGCCAATTGAAG-1      1
chr1    10006   10386   TCGTTATTCATCCTGC-1      1
chr1    10006   10386   TGTTCCTCAATAGCAA-1      1
chr1    10006   10443   CGCCTGTGTTGGCGTG-1      1
chr1    10007   10013   CGTATTGCACTAGGTC-1      1
chr1    10007   10018   AGCTTTAAGGTCGAGG-1      1
chr1    10007   10018   GCCCATAAGCACGATT-1      1
chr1    10007   10019   CGTCATTGTGTCACGG-1      1
```
Cell-Type File: 
```bash
#csv file with barcode and cell type. Here, cancer cells are labeled as "Cancer" but the cancer label can be set using the cancerLabel parameter in sort_fragments()

$ head celltypes.csv
GGTTTAATCTAAGGTC-1,Normal
GGTTTCCTCAATCTCT-1,Cancer
GGTTTCCTCATAAGCC-1,Cancer
GGTTTCCTCTAATTGG-1,Cancer
GGTTTGTAGCGATACT-1,Cancer
GGTTTGTAGGGATGAC-1,Normal
GTAAAGCCACTTACAG-1,Cancer
GTAAGCGCAAGCGATG-1,Cancer
GTAAGCGCACCGGCTA-1,Cancer
GTAAGCGCATAAGCAA-1,Cancer
GTAAGCTTCGAAGTAG-1,Normal
```



### Parameters

| Parameter | Description | Default |
|------------|-------------|----------|
| `window_size` | Number of cancer reads per window. Values between 2 and 100 work best.  | `10` |
| `step_size` | Step size between consecutive windows | `1` |
| `z_cutoff` | Z-score cutoff for significance. Values between 4 and 5 work best. | `4.5` |
| `chromosomes` | List of chromosomes to analyze | All human chromosomes |


### Output Format

The resulting file (`deletions.tsv`) contains:
| Column | Description |
|---------|-------------|
| `chrom` | Chromosome name |
| `start` | Start position of deletion region |
| `end` | End position of deletion region |
| `mean_pval` | Mean p-value of merged windows within region |
| `min_pval` | Minimum p-value of merged windows within region |
| `max_pval` | Maximum p-value of merged windows within region |
| `C_counts` | Read count in cancer cells|
| `NC_counts` | Read count in non-cancer cells |
