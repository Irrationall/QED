from dataclasses import dataclass, asdict, field
import pandas as pd
import pkg_resources
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
    



def merge_df(geneset_list: List[geneset], term_suffix: bool = True) :

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
        
        df['Term'] = df['Term'].str.strip()

    return df 




# DB list

@dataclass
class EnrichR_DB:

    db: list = field(default_factory=list)

    def add_info(self, info):
        self.db.append(info)

    def __repr__(self) :
        
        return f"ENRICHR DATABASE\n[number of databases: {len(self.db)}]"
    
    @property
    def categories(self):
        return list(set(info.category for info in self.db))

    def search_DB(self, 
                  category: List[str] = None, 
                  is_latest: bool = True,
                  by_name: str = None):
        
        if is_latest == True:
            db_list = [db for db in self.db if db.category in category and db.is_latest == True]
            
            if by_name is not None:
                db_list = [db for db in db_list if by_name.lower() in db.name.lower()]

        else:
            db_list = [db for db in self.db if db.category in category]
            
            if by_name is not None:
                db_list = [db for db in db_list if by_name.lower() in db.name.lower()]

        return db_list 




@dataclass
class EnrichR_DB_info :
    
    name: str
    category: str
    description: str
    number_of_terms: int
    gene_coverage: int
    genes_per_term: int
    is_latest: bool 
    
    db: EnrichR_DB = field(default_factory=EnrichR_DB)

    def __post_init__(self):
        self.db.add_info(self)

    def __repr__(self) :
        
        return f"name: {self.name} / category: {self.category}\n"


ENRICHR_DB = EnrichR_DB()

# file contains EnrichR database information
file_path = pkg_resources.resource_filename('qed.data', 'DB/EnrichR_DB_list.csv')
DB_DF = pd.read_csv(file_path)

for _, row in DB_DF.iterrows():
    EnrichR_DB_info(
        name = row['DB'],
        category = row['Category'],
        description =  row['Description'],
        number_of_terms = row['Terms'],
        gene_coverage=row['Gene coverage'],
        genes_per_term=row['Genes per term'],
        is_latest=row['is_latest'],
        db=ENRICHR_DB
    )

