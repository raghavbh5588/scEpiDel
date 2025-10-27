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