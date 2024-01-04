import matplotlib.pyplot as plt
from typing import Dict, Tuple
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
            cmap: str = 'Blues',
            xlabel_rotation: int = 90,
            figsize: Tuple = (5,10),
            vmin: float = None,
            vmax: float = None,
            cbar_kws: Dict = None) :
        
        """Draw heatmap with merged DataFrame

        Args:
            df (pd.DataFrame): Dataframe that contains EnrichR Query result
            
            n (int): Number of terms for each 'group_by'
            
            group_by (str): A column from Dataframe for x axis
            
            order_by (str, optional): Determine criteria when select terms. 
                                      Defaults to "Adjusted p-value".
                                      
            allow_duplicate (bool, optional): Allow or not duplicates in terms. 
                                              Defaults to False.
                                              
            method (str, optional): Method of scipy.cluster.hierarchy.linkage. Defaults to 'average'.
            
            metric (str, optional): Metric of scipy.cluster.hierarchy.linkage. Defaults to 'euclidean'.
            
            cmap (str, optional): colormap of heatmap. Passed to matplotlib.pyplot.pcolormesh. 
                                  Defaults to 'Blues'.
                                  
            xlabel_ratation (int, optional): Angle that x axis labels will be rotated
            
            figsize (Tuple, optional): Determine figure size. Passed to matplotlib.pyplot.subplots()
            
            vmin, vmax (float, optional): Determine data range that colormap covers
                                          Defaults to None  
                               
            cbar_kws (Dict, optional): colorbar arguments. Passed to matplotlib.pyplot.colorbar                      

        Returns:
            matplotlib.axes class: You can manage plot with matplotlib package
        """

        if order_by not in ['Adjusted p-value', 'P-value', 'Odds ratio', 'Combined score'] :
             raise ValueError("'Order_by' should be one of 'Adjusted p-value', 'P-value', 'Odds ratio', and 'Combined score'.")


        fig, ax = plt.subplots(figsize=figsize)
        
        subset_df = select_top_n(df, 
                                 n, 
                                 group_by, 
                                 order_by, 
                                 allow_duplicate)
        
        if order_by in ['Adjusted p-value', 'P-value'] :

            pivot_df = subset_df.pivot(index = "Term", columns = group_by, values = order_by).fillna(1).rename_axis(None, axis=1)
            pivot_df = - np.log10(pivot_df)

        elif order_by in ['Odds ratio', 'Combined score'] :
            pivot_df = subset_df.pivot(index = "Term", columns = group_by, values = order_by).rename_axis(None, axis=1)

        row_linkage = hierarchy.linkage(pivot_df.values, method=method, metric=metric)
        row_order = hierarchy.leaves_list(row_linkage)

        col_linkage = hierarchy.linkage(pivot_df.values.T, method=method, metric=metric)
        col_order = hierarchy.leaves_list(col_linkage)

        clustered_data = pivot_df.iloc[row_order, col_order]
        clustered_data

        quadmesh = ax.pcolormesh(clustered_data, 
                                 cmap=cmap, 
                                 vmin=vmin, 
                                 vmax=vmax)

        plt.xticks(np.arange(len(clustered_data.columns)) + 0.5, clustered_data.columns, rotation=xlabel_rotation, ha='center')
        plt.yticks(np.arange(len(clustered_data.index)) + 0.5, clustered_data.index)
        
        if cbar_kws is None :
            cbar = fig.colorbar(quadmesh, ax=ax, label=f'-log10{order_by}', shrink = 0.5)
        else :
            cbar = fig.colorbar(quadmesh, ax=ax, **cbar_kws)
        
        plt.gca().invert_yaxis()

        return fig, ax, cbar