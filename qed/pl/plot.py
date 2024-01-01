import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
import numpy as np
import pandas as pd
from .organize import select_top_n




def heatmap(df: pd.DataFrame,
            n: int,
            group_by: str,
            order_by: str = "Adjusted p-value",
            allow_duplicate: bool = False,
            method: str = 'average',
            metric: str = 'euclidean',
            cmap: str = 'Blues') :

    subset_df = select_top_n(df, n, group_by, order_by, allow_duplicate)

    pivot_df = subset_df.pivot(index = "Term", columns = group_by, values = order_by).fillna(1).rename_axis(None, axis=1)
    pivot_df = - np.log10(pivot_df)

    row_linkage = hierarchy.linkage(pivot_df.values, method=method, metric=metric)
    row_order = hierarchy.leaves_list(row_linkage)

    col_linkage = hierarchy.linkage(pivot_df.values.T, method=method, metric=metric)
    col_order = hierarchy.leaves_list(col_linkage)

    clustered_data = pivot_df.iloc[row_order, col_order]
    clustered_data

    ax = plt.pcolormesh(clustered_data, cmap=cmap)

    plt.xticks(np.arange(len(clustered_data.columns)) + 0.5, clustered_data.columns, rotation=90, ha='center')
    plt.yticks(np.arange(len(clustered_data.index)) + 0.5, clustered_data.index)

    plt.colorbar(label='-log10(adjusted p-value)', shrink = 0.5)
    plt.gca().invert_yaxis()

    return ax