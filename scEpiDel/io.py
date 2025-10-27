import pandas as pd
from tqdm import tqdm
import numpy as np
import os

def rmoverlap(df, chroms):
    keep_indices = []
    for chrom in tqdm(chroms, desc="Removing Overlapping Reads..."):
        sub_df = df[df[0] == chrom]
        starts = sub_df[1].to_numpy()
        ends = sub_df[2].to_numpy()
        keep_mask = np.zeros(len(sub_df), dtype=bool)
        prev_end = -1
        for i in range(len(sub_df)):
            if starts[i] >= prev_end:
                keep_mask[i] = True
                prev_end = ends[i]
        keep_indices.extend(sub_df.index[keep_mask])
    return df.loc[keep_indices].reset_index(drop=True)


def sort_fragments(fragmentFile, cellTypeFile, cancerLabel, outputDir, chromosomes):
    # Ensure output directory exists
    os.makedirs(outputDir, exist_ok=True)

    cellTypes = pd.read_csv(cellTypeFile, header=None)
    barcodes_cancer = set(cellTypes.loc[cellTypes.iloc[:, 1] == cancerLabel, cellTypes.columns[0]])
    barcodes_normal = set(cellTypes.loc[cellTypes.iloc[:, 1] != cancerLabel, cellTypes.columns[0]])

    if len(barcodes_cancer) == 0:
        raise ValueError("No cancer barcodes found with the specified cancerLabel.")
    if len(barcodes_normal) == 0:
        raise ValueError("No normal barcodes found in the cell type file.")
    
    cancer_frags = []
    normal_frags = []
    
    print("Reading and sorting fragment file...")
    with open(fragmentFile, "r") as fragments:
        for line in fragments:
            parts = line.split('\t')
            barcode = parts[3].strip()
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            
            if barcode in barcodes_cancer:
                cancer_frags.append((chrom, start, end, barcode))
            elif barcode in barcodes_normal:
                normal_frags.append((chrom, start, end, barcode))
    
    if len(cancer_frags) == 0:
        raise ValueError("No cancer fragments found for the specified cancerLabel.")
    if len(normal_frags) == 0:
        raise ValueError("No normal fragments found in the fragment file.")
    
    print("Removing overlaps from cancer fragments...")
    cancer_df = pd.DataFrame(cancer_frags, columns=[0, 1, 2, 3])
    cancer_df_filtered = rmoverlap(cancer_df, chromosomes)
    
    print("Removing overlaps from normal fragments...")
    normal_df = pd.DataFrame(normal_frags, columns=[0, 1, 2, 3])
    normal_df_filtered = rmoverlap(normal_df, chromosomes)

    print("Writing filtered fragments...")

    cancerFile = os.path.join(outputDir, "fragments_cancer_filtered.bed")
    normalFile = os.path.join(outputDir, "fragments_normal_filtered.bed")

    cancer_df_filtered.to_csv(cancerFile, sep='\t', header=False, index=False)
    normal_df_filtered.to_csv(normalFile, sep='\t', header=False, index=False)
    
    return cancerFile, normalFile
