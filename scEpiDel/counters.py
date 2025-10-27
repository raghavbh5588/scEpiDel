import numpy as np
import pandas as pd
from tqdm import tqdm

def counter(read_array, start, end):
    left_index = np.searchsorted(read_array, start, side='right')
    right_index = np.searchsorted(read_array, end, side='left')
    count = right_index - left_index
    return count

def count_reads(C_df, NC_df, chroms, window_size, step_size):
    counts_results = []
    for chrom in tqdm(chroms, desc = "Scanning genome..."):
        C_sub_df = C_df[C_df[0] == chrom]
        C_starts = C_sub_df[1].to_numpy()
        C_ends = C_sub_df[2].to_numpy()
        
        NC_sub_df = NC_df[NC_df[0] == chrom]
        NC_mid = ((NC_sub_df[2] + NC_sub_df[1]) / 2).to_numpy()
        
        n_windows = (len(C_sub_df) - window_size) // step_size + 1
        
        for i in range(n_windows):
            s_index = i * step_size
            e_index = s_index + window_size - 1
            
            NC_count = counter(NC_mid, C_starts[s_index], C_ends[e_index])
            
            counts_results.append([chrom, C_starts[s_index], C_ends[e_index], NC_count])
    
    counts_df = pd.DataFrame(counts_results, columns=['chrom', 'start', 'end', 'NC_count'])
    return counts_df