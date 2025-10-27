import pandas as pd
import numpy as np
from scipy.stats import norm
from .counters import count_reads
from .merge import merge_windows

def scEpiDel(cancer_file, normal_file, chromosomes, step_size=1, window_size=10, z_cutoff=4.5):
    print("Loading fragments...")
    C_df = pd.read_csv(cancer_file, sep="\t", header=None)
    NC_df = pd.read_csv(normal_file, sep="\t", header=None)

    print("Counting reads...")
    count_df = count_reads(C_df, NC_df, chromosomes, window_size, step_size)
    x = np.log(count_df["NC_count"] + 1)
    count_df["NC_count_Z"] = (x - np.mean(x)) / np.std(x)
    count_df["p_value"] = 1 - norm.cdf(count_df["NC_count_Z"])
    filtered = count_df[count_df["NC_count_Z"] > z_cutoff]

    if len(filtered) == 0:
        print("No significant windows found.")
        return pd.DataFrame(columns=['chrom', 'start', 'end', 'mean_pval', 'C_counts', 'NC_counts'])

    print(f"Found {len(filtered)} significant windows.")
    merged = merge_windows(filtered, chromosomes, C_df, NC_df)
    return merged
