import pandas as pd




def select_top_n(df: pd.DataFrame, 
                 n: int,
                 group_by: str,
                 order_by: str = 'Adjusted p-value', 
                 allow_duplicate: bool = False):
    
    # Check if 'Term' is in dataframe

    if 'Term' not in df.columns :
        raise ValueError("The 'Term' column is not present in the input DataFrame.")
    
    # Chcek 'order_by' is valid one.

    if order_by in ['Adjusted p-value', 'P-value'] :
        ascending = True
    elif order_by in ['Odds ratio', 'Combined score'] :
        ascending = False
    else :
        raise ValueError("'Order_by' should be one of 'Adjusted p-value', 'P-value', 'Odds ratio', and 'Combined score'.")
    
    df_copy = df.copy(deep=True)
    df_copy = df_copy.sort_values([group_by, order_by], ascending=ascending)
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
            
            df_top_n = df_top_n.sort_values(order_by, ascending=ascending)
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

    # Return 'Overlapping genes' to list    

    if 'Overlapping genes' in df_top_n.columns :
        df_top_n['Overlapping genes'] = df_top_n['Overlapping genes'].apply(list)

    # Subset original df using 'Term' and 'Database'

    df_subset = pd.DataFrame()

    for _, row in df_top_n.iterrows():

        # Extract 'Term' and 'Database' 

        term = row['Term']
        database = row['Database']

        df_row_subset = df[(df['Term'] == term) & (df['Database'] == database)]
        
        df_subset = pd.concat([df_subset, df_row_subset])

    df_subset = df_subset.reset_index(drop=True)

    return df_subset