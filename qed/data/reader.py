from .structure import geneset
from typing import List, Optional, TYPE_CHECKING
import pandas as pd
import time

if TYPE_CHECKING:
    from anndata import AnnData




# Wrapper for calculate time
def calc_time(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"Execution time for {func.__name__}: {end_time - start_time} seconds")
        return result
    return wrapper




@calc_time
def readtxt(file_path: str, sep: str, format: Optional[str] = 'rowside') -> List[geneset] :
    
    setlist = []

    if format == 'rowside' :
        with open(file_path, 'r') as f :
            lines = [line.rstrip() for line in f]

        for line in lines :
            name = line.split(sep)[0]
            genes = list(filter(None, line.split(sep)[1:]))
            Gene_Set = geneset(name, genes)
            setlist.append(Gene_Set)
    
    elif format == 'colside' :
        df = pd.read_csv(file_path, sep=sep)
        colnames = df.columns

        for col in colnames :
             name = col
             genes = df[col].dropna().tolist()
             Gene_Set = geneset(name, genes)
             setlist.append(Gene_Set)

    else :
        raise ValueError("Invalid 'format' parameter. Supported values: 'rowside', 'colside'")

    return setlist


@calc_time
def readseurat(file_path: str, sep: str, 
               adj_pval_cutoff: float = 0.05,
               lfc_cutoff: float = 0.5,
               pct_cutoff: float = 0.05, 
               *args, **kwargs) -> List[geneset] :
    
    setlist = []
    
    df = pd.read_csv(file_path, sep=sep, *args, **kwargs)
    df = df[(df['p_val_adj'] <= adj_pval_cutoff) & (df['avg_log2FC'] >= lfc_cutoff) & (df['pct.1'] >= pct_cutoff)]
    
    tmp_group = df.groupby('cluster')['gene']
    
    for name, group in tmp_group:
        genes = group.tolist()
        gs = geneset(name=name, genes=genes)
        setlist.append(gs)
    
    return setlist




@calc_time
def readscanpy(adata: 'AnnData',
               adj_pval_cutoff: float = 0.05,
               lfc_cutoff: float = 0.5,
               select_top_n: Optional[int] = None,
               select_order: Optional[str] = 'scores') -> List[geneset] :
    
    
    if 'rank_genes_groups' not in adata.uns:
        raise ValueError("Run scanpy.tl.rank_genes_groups first.")
    
    setlist = []

    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names

    new_df = pd.DataFrame()

    for group in groups:
        
        temp_df = pd.DataFrame()
        temp_df['cluster_names'] = pd.Series([group]*len(result['names'][group]))
        temp_df['names'] = result['names'][group]
        temp_df['scores'] = result['scores'][group]
        temp_df['pvals'] = result['pvals'][group]
        temp_df['pvals_adj'] = result['pvals_adj'][group]
        temp_df['logfoldchanges'] = result['logfoldchanges'][group]
        
        new_df = pd.concat([new_df, temp_df])

    new_df.reset_index(drop=True, inplace=True)

    new_df = new_df[(new_df['logfoldchanges'] > lfc_cutoff) &
                    (new_df['pvals_adj'] < adj_pval_cutoff)
                   ]
    
    if select_top_n is not None :
        
        top_df = new_df.groupby('cluster_names').apply(lambda x: x.nlargest(select_top_n, select_order)).reset_index(drop=True)
        tmp_group = top_df.groupby('cluster_names')['names']

        for name, group in tmp_group:
            
            genes = group.tolist()
            gs = geneset(name=name, genes=genes)
            setlist.append(gs)

    
    else :
        
        tmp_group = new_df.groupby('cluster_names')['names']

        for name, group in tmp_group:
            
            genes = group.tolist()
            gs = geneset(name=name, genes=genes)
            setlist.append(gs)
     
    return setlist