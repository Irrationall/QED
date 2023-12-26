from dataclasses import dataclass, asdict, field
import pandas as pd
from typing import List, Dict


@dataclass
class geneset:
    name: str
    genes: List[str]
    annotation: str = None
    GO: List[pd.DataFrame] = field(default_factory=list)
    params: Dict = field(default_factory=dict)
    
    
    def __repr__(self) :
        return f"geneset object [name: {self.name}, number of genes: {len(self.genes)}, number of GO_dataframes: {len(self.GO)}]"
    
    def resetGO(self) :
        self.GO = []
        return self
    



def merge_df(geneset_list, term_suffix=True) :

    merged_dfs = []

    for geneset in geneset_list:
        merged_df = pd.concat(geneset.GO, axis=0, ignore_index=True)
        merged_dfs.append(merged_df)

    df = pd.concat(merged_dfs, axis=0)

    if term_suffix == False :
        # Remove (GO:000) patterns at the end of the 'Term' column for Gene Ontology
        df['Term'] = df['Term'].str.replace(r'\(GO:\d+\)$', '', regex=True)

        # Remove R-HSA-000 patterns at the end of the 'Term' column for Reactome
        df['Term'] = df['Term'].str.replace(r'\bR-HSA-\d+\b', '', regex=True)

        # Remove WP000 patterns at the end of the 'Term' column for WikiPathway
        df['Term'] = df['Term'].str.replace(r'\bWP\d+\b', '', regex=True)

        # Remove MP:0000 at the end of the 'Term' column for MGI Mammalian Phenotype Level 4
        df['Term'] = df['Term'].str.replace(r'MP:\d+$', '', regex=True)

        # Remove (HP:0000) at the end of the 'Term' column for Human phenotype Ontology
        df['Term'] = df['Term'].str.replace(r'\(HP:\d+\)$', '', regex=True)

        # Remove GSE000 at the end of the 'Term' column for several database
        df['Term'] = df['Term'].str.replace(r'GSE\d+$', '', regex=True)

        # Remove GSD000 at the end of the 'Term' column for several database
        df['Term'] = df['Term'].str.replace(r'GSD\d+$', '', regex=True)

        # Remove CL:000 at the end of the 'Term' column for Tabula Muris
        df['Term'] = df['Term'].str.replace(r'CL:\d+$', '', regex=True)

        # Remove CL000 at the end of the 'Term' column for Azimuth Celltype 2021
        df['Term'] = df['Term'].str.replace(r'CL\d+$', '', regex=True)

    return df 




