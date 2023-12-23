from typing import List, Dict, Any
from .structure import geneset
import concurrent.futures
from concurrent.futures import as_completed
from threading import Lock
from tqdm.notebook import tqdm
import pandas as pd
import requests
import json




def upload_genes(genes: List[str]) :

    """ Upload gene sets to EnrichR website """

    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
    genes_str = '\n'.join(genes)
    description = 'Example gene list'
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }

    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')

    data = json.loads(response.text)

    return data


def to_dataframe(request_res: Dict, database: str, annot_colname: str, annot: Any) -> pd.DataFrame:

    """ Convert response of request to pandas dataframe
        * request result = JSON format
        * returns pandas.DataFrame """

    colnames = ["Rank", "Term", "P-value", "Odds ratio", "Combined score", "Overlapping genes", "Adjusted p-value", "Old p-value", "Old adjusted p-value"]
    df = pd.DataFrame(request_res[database], columns=colnames)
    df['Database'] = database
    df[annot_colname] = annot

    return df


def get_enrichment_data(genes: List, database: str, annot_colname: str, annot: Any) :

    """ Get response from EnrichR website using querying gene set"""
    
    data = upload_genes(genes)

    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'
    user_list_id = data['userListId']
    gene_set_library = database

    url = ENRICHR_URL + query_string % (user_list_id, gene_set_library)
    response = requests.get(url, stream=True)
    if not response.ok:
        print(f"Status code: {response.status_code}")
        print(f"Response content: {response.content}")
        raise Exception('Error fetching enrichment results')
    
    res = json.loads(response.text)

    df = to_dataframe(res, database, annot_colname, annot)
    
    return df


def _get_multiple_enrichment_data(geneset: geneset, database: str, annot_colname: str, annot: Any = None):
    
    try:
        annot = annot if annot is not None else geneset.name
        result_df = get_enrichment_data(geneset.genes, database, annot_colname, annot)
        geneset.GO.append(result_df)
    
    except Exception as e:
        print(f"Error processing {geneset.name} for {database}: {str(e)}")

    return geneset


def get_enrichment_dataframes_spare(geneset_list :List[geneset], dblist: List, annot_colname: str, annot: Any = None):
    print_lock = Lock()

    with concurrent.futures.ThreadPoolExecutor() as executor:
        for geneset in geneset_list:
            futures = {executor.submit(_get_multiple_enrichment_data, geneset, database, annot_colname, annot): database for database in dblist}
            concurrent.futures.wait(futures)
            
            with print_lock:
                print(f'{geneset.name} completed!')
    
    with print_lock :    
        print(f'All process completed!')

    return geneset_list



def get_enrichment_dataframes(geneset_list: List[geneset], dblist: List, annot_colname: str, annot: Any =None, n_jobs: int = None):

    with concurrent.futures.ThreadPoolExecutor(max_workers = n_jobs) as executor:

        futures = {executor.submit(_get_multiple_enrichment_data, geneset, database, annot_colname, annot): geneset for geneset in geneset_list for database in dblist}
        
        pbar = tqdm(total=len(geneset_list), desc='Processing genesets')
        
        for future in as_completed(futures):
            geneset = futures[future]
            pbar.update(1)
        
    pbar.close()

    return geneset_list








