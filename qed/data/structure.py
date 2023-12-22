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


@dataclass
class enrichrdb:
    name: str
    short: str


def merging_df(geneset_list) :

    merged_dfs = []

    for geneset in geneset_list:
        merged_df = pd.concat(geneset.GO, axis=0, ignore_index=True)
        merged_dfs.append(merged_df)

    df = pd.concat(merged_dfs, axis=0)

    return df 


