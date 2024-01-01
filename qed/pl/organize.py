import pandas as pd




def select_top_n(df: pd.DataFrame, 
                 n: int,
                 group_by: str,
                 order_by: str = 'Adjusted p-value', 
                 allow_duplicate: bool = False):
    
    if 'Term' not in df.columns :
        raise ValueError("The 'Term' column is not present in the input DataFrame.")
    
    df_copy = df.copy(deep=True)
    df_copy = df_copy.sort_values([group_by, order_by], ascending=True)
    df_top_n = df_copy.groupby(group_by).head(n)
    

    if not allow_duplicate :

        if 'Overlapping genes' in df_copy.columns:

            df_copy = df_copy.copy()
            df_copy['Overlapping genes'] = df_copy['Overlapping genes'].apply(tuple)
            df_top_n = df_top_n.copy()
            df_top_n['Overlapping genes'] = df_top_n['Overlapping genes'].apply(tuple)

        while df_top_n['Term'].duplicated().any():
            
            # Delete rows in df_copy that presents in df_top_n
            
            duplicates = df_top_n['Term'][df_top_n['Term'].duplicated()].unique()
            print(f'***Duplicated Terms***\n{duplicates}')

            df_copy = pd.merge(df_copy, df_top_n, indicator=True, how='outer').query('_merge=="left_only"').drop('_merge', axis=1)

            # Find and remove duplicates
            
            df_top_n = df_top_n.sort_values(order_by)
            mask = df_top_n.duplicated(subset='Term', keep='first')
            
            deleted_rows = {}
            for celltype in df_top_n.loc[mask, group_by].unique():
                deleted_rows[celltype] = (mask & (df_top_n[group_by] == celltype)).sum()
                
            df_top_n = df_top_n[~mask]

            # Refill df_top_n
            
            for celltype, count in deleted_rows.items():
                next_terms = df_copy[df_copy[group_by] == celltype].nsmallest(count, order_by)
                df_top_n = pd.concat([df_top_n, next_terms])
                
    df_top_n = df_top_n.sort_values([group_by, order_by])

    if 'Overlapping genes' in df_top_n.columns :
        df_top_n['Overlapping genes'] = df_top_n['Overlapping genes'].apply(list)

    df_top_n = df_top_n.reset_index(drop=True)

    return df_top_n