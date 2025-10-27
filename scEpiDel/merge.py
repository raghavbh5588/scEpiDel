import numpy as np
import pandas as pd
from tqdm import tqdm

def merge_windows(df, chroms, C_reads, NC_reads):
    all_merged_dfs = []

    for chrom in tqdm(chroms, desc="Merging windows..."):
        sub_df = df[df.iloc[:, 0] == chrom].sort_values(by=df.columns[1])
        starts = sub_df.iloc[:, 1].to_numpy()
        ends = sub_df.iloc[:, 2].to_numpy()
        pvals = sub_df["p_value"].to_numpy()

        if len(starts) == 0:
            continue

        C_sub_df = C_reads[C_reads.iloc[:, 0] == chrom]
        NC_sub_df = NC_reads[NC_reads.iloc[:, 0] == chrom]

        C_mid = np.sort((C_sub_df.iloc[:, 1] + C_sub_df.iloc[:, 2]) / 2)
        NC_mid = np.sort((NC_sub_df.iloc[:, 1] + NC_sub_df.iloc[:, 2]) / 2)

        merged_windows_chrom = []
        win_start = starts[0]
        win_end = ends[0]
        current_pvals = [pvals[0]]

        for i in range(1, len(starts)):
            if starts[i] <= win_end:
                win_end = max(win_end, ends[i])
                current_pvals.append(pvals[i])
            else:
                merged_windows_chrom.append([
                    chrom,
                    win_start,
                    win_end,
                    np.mean(current_pvals),
                    np.min(current_pvals),
                    np.max(current_pvals)
                ])
                win_start = starts[i]
                win_end = ends[i]
                current_pvals = [pvals[i]]

        # append the last region
        merged_windows_chrom.append([
            chrom,
            win_start,
            win_end,
            np.mean(current_pvals),
            np.min(current_pvals),
            np.max(current_pvals)
        ])

        merged_df = pd.DataFrame(
            merged_windows_chrom,
            columns=['chrom', 'start', 'end', 'mean_pval', 'min_pval', 'max_pval']
        )

        # Count reads
        left_indices_C = np.searchsorted(C_mid, merged_df["start"].to_numpy(), side='right')
        right_indices_C = np.searchsorted(C_mid, merged_df["end"].to_numpy(), side='left')
        merged_df["C_counts"] = right_indices_C - left_indices_C

        left_indices_NC = np.searchsorted(NC_mid, merged_df["start"].to_numpy(), side='right')
        right_indices_NC = np.searchsorted(NC_mid, merged_df["end"].to_numpy(), side='left')
        merged_df["NC_counts"] = right_indices_NC - left_indices_NC

        all_merged_dfs.append(merged_df)

    final_df = pd.concat(all_merged_dfs, ignore_index=True)
    return final_df
