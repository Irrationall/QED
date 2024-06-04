from typing import List, Dict, Any
from .structure import geneset
import concurrent.futures
from concurrent.futures import as_completed
from threading import Lock
from tqdm import tqdm
import pandas as pd
import requests
import json
import copy




# Basic Gene Ontoloy Analysis

def upload_genes(genes: List[str]) :

    """ Upload gene sets to EnrichR website

         Args
            genes (List) : A list of genes to query to enrichR

    """

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




def to_dataframe(request_res: Dict, 
                 database: str, 
                 annot_colname: str, 
                 annot: Any) -> pd.DataFrame:

    """ Convert response of request to pandas dataframe
        * request result = JSON format
        * returns pandas.DataFrame 
    """

    colnames = ["Rank", "Term", "P-value", "Odds ratio", "Combined score", "Overlapping genes", "Adjusted p-value", "Old p-value", "Old adjusted p-value"]
    df = pd.DataFrame(request_res[database], columns=colnames)
    df['Database'] = database
    df[annot_colname] = annot

    return df




def get_enrichment_data(genes: List, 
                        database: str, 
                        annot_colname: str, 
                        annot: Any) :

    """ Get response from EnrichR website using querying gene set
    
         Args
            genes (List): A list of genes to query for enrichR.

            database (str): Name of database to query for enrichR.

            annot_colname (str): A column name for each gene sets.

            annot (Any): value for column 'annot_colname'.

    """
    
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




def _get_multiple_enrichment_data(geneset: geneset, 
                                  database: str, 
                                  annot_colname: str, 
                                  annot: Any = None):
    
    try:
        annot = annot if annot is not None else geneset.name
        result_df = get_enrichment_data(geneset.genes, database, annot_colname, annot)
        geneset.GO.append(result_df)

        geneset.params['annot_colname'] = annot_colname
        geneset.params['annot'] = annot
    
    except Exception as e:
        
        print(f"Error processing {geneset.name} for {database}: {str(e)}")

    return geneset




def get_enrichment_dataframes_spare(geneset_list :List[geneset], 
                                    dblist: List, annot_colname: str, 
                                    annot: Any = None):
    
    print_lock = Lock()

    with concurrent.futures.ThreadPoolExecutor() as executor:
        for geneset in geneset_list:
            futures = {executor.submit(_get_multiple_enrichment_data, 
                                       geneset, 
                                       database, 
                                       annot_colname, 
                                       annot): database for database in dblist}
            concurrent.futures.wait(futures)
            
            with print_lock:
                print(f'{geneset.name} completed!')
    
    with print_lock :    
        print(f'All process completed!')

    return geneset_list




def get_enrichment_dataframes(geneset_list: List[geneset], 
                              dblist: List, 
                              annot_colname: str, 
                              annot: Any = None, 
                              n_jobs: int = None,
                              handle_error: bool = False,
                              max_iter: int = 10):
    
    _geneset_list = copy.deepcopy(geneset_list)

    for geneset in _geneset_list:
        geneset.params['database'] = dblist
    
    """
    Retrieves enrichment dataframes for a list of genesets from multiple databases.

    Args:
        geneset_list (List[geneset]): A list of genesets for enrichment analysis.

        dblist (List): A list of databases to perform enrichment analysis on.

        annot_colname (str): The column name in the database for annotation.

        annot (Any, optional): Additional annotation information. Defaults to None.

        max_iter (int, optional): The maximum number of iterations to run. Defaults to 10.
        
        n_jobs (int, optional): The number of parallel jobs to run. Defaults to None.

    Returns:
        List[geneset]: The input geneset list.

    """
    with concurrent.futures.ThreadPoolExecutor(max_workers=n_jobs) as executor:
        futures = {executor.submit(_get_multiple_enrichment_data, 
                                   geneset, 
                                   database, 
                                   annot_colname, 
                                   annot): geneset for geneset in _geneset_list for database in dblist}
        
        pbar = tqdm(total=len(geneset_list), desc='Processing genesets')
        
        for future in concurrent.futures.as_completed(futures):

            geneset = futures[future]
            pbar.update(1)
        
    pbar.close()


    # Re-request for geneset that has error
    if handle_error :

        print('Trying to re-request...')

        geneset_list_2 = copy.deepcopy(_geneset_list)

        for geneset in geneset_list_2:
            
            len_success = len(geneset.GO)
            db_successs = [df['Database'].values[0] for df in geneset.GO]
            retry_db = [db for db in geneset.params['database'] if db not in db_successs]

            if len_success < len(geneset.params['database']) :
                geneset.params['re_request'] = retry_db
            else :
                geneset.params['re_request'] = []

        if all([len(geneset.params['re_request']) == 0 for geneset in geneset_list_2]):
            print("Nothing to re-request")

        max_iter = max_iter
        iter_count = 0

        # Copilot version. Have to handle it.
        while any([len(geneset.params['re_request']) > 0 for geneset in geneset_list_2]) and iter_count < max_iter :

            for geneset in geneset_list_2:
                
                for database in geneset.params['re_request']:
                    
                    _get_multiple_enrichment_data(geneset, database, annot_colname, annot)

                    print(f'{geneset.name} with {database} completed!')

            for geneset in geneset_list_2:
                
                len_success = len(geneset.GO)
                db_successs = [df['Database'].values[0] for df in geneset.GO]
                retry_db = [db for db in geneset.params['re_request'] if db not in db_successs]

                if len_success < len(geneset.params['database']) :
                    geneset.params['re_request'] = retry_db
                else :
                    geneset.params['re_request'] = []

            iter_count += 1

    return geneset_list_2 if handle_error else _geneset_list




# Gene Ontology Analysis with background genes

def upload_genes_with_background (genes: List[str]) :
    
    base_url = "https://maayanlab.cloud/speedrichr"

    genes = genes.copy()

    description = "sample gene set with background"

    res = requests.post(
        base_url+'/api/addList',
        files=dict(
        list=(None, '\n'.join(genes)),
        description=(None, description),
        )
    )
    if res.ok:
        userlist_response = res.json()
        
    return userlist_response['userListId']




def upload_background_genes(genes: List[str]) :

    base_url = "https://maayanlab.cloud/speedrichr"

    bggenes = genes.copy()

    res = requests.post(
        base_url+'/api/addbackground',
        data=dict(background='\n'.join(bggenes)),
    )

    if res.ok:
        background_response = res.json()

    return background_response['backgroundid']




def get_enrichment_data_with_background(query_genes: List[str],
                                        background_genes: List[str], 
                                        database: str, 
                                        annot_colname: str, 
                                        annot: Any
                                        ) :

    """ Get response from EnrichR website using querying and background gene set
    
        Args:
            query_genes (List): List of genes want to query

            background_genes (List) : List of gene want to be set to background

            database (str) : Name of database for backgroundType  
    
    """
    
    UID = upload_genes_with_background(query_genes)
    BGID = upload_background_genes(background_genes)

    base_url = "https://maayanlab.cloud/speedrichr"

    res = requests.post(
        base_url+'/api/backgroundenrich',
        data = dict(
            userListId = UID,
            backgroundid = BGID,
            backgroundType = database,
        )
    )

    if res.ok:
        results = res.json()

    df = to_dataframe(results, database, annot_colname, annot)
    
    return df




def _get_multiple_enrichment_data_with_background(
    geneset: geneset,
    background_geneset: List,
    database: str,
    annot_colname: str,
    annot: Any = None
):
    """
    Get multiple enrichment data with background geneset.

    Args:
        geneset (geneset): The geneset object containing the genes of interest.

        background_geneset (List): The list of background genes.

        database (str): The name of the database to query.

        annot_colname (str): The name of the column containing annotations in the database.

        annot (Any, optional): The annotation to use. Defaults to None.

    Returns:
        geneset: The updated geneset object with enrichment data appended to the GO attribute.
    """
    try:
        annot = annot if annot is not None else geneset.name

        result_df = get_enrichment_data_with_background(
            geneset.genes,
            background_geneset,
            database,
            annot_colname,
            annot
        )
        geneset.GO.append(result_df)

    except Exception as e:
        print(f"Error processing {geneset.name} for {database}: {str(e)}")

    return geneset




def get_enrichment_dataframes_with_background_spare(geneset_list :List[geneset], 
                                    dblist: List, annot_colname: str, 
                                    annot: Any = None):
    
    print_lock = Lock()

    with concurrent.futures.ThreadPoolExecutor() as executor:
        for geneset in geneset_list:
            futures = {executor.submit(_get_multiple_enrichment_data_with_background, 
                                       geneset, 
                                       database, 
                                       annot_colname, 
                                       annot): database for database in dblist}
            concurrent.futures.wait(futures)
            
            with print_lock:
                print(f'{geneset.name} completed!')
    
    with print_lock :    
        print(f'All process completed!')

    return geneset_list




def get_enrichment_dataframes_with_background(geneset_list: List[geneset],
                              background_geneset: List, 
                              dblist: List, 
                              annot_colname: str, 
                              annot: Any =None, 
                              n_jobs: int = None):
    """
    Get enrichment dataframes with background genes.

    Args:
        geneset_list (List[geneset]): A list of genesets.
        
        background_geneset (List): A list of background genes.
        
        dblist (List): A list of databases.
        
        annot_colname (str): The name of the annotation column.
        
        annot (Any, optional): Additional annotation. Defaults to None.
        
        n_jobs (int, optional): The number of parallel jobs to run. Defaults to None.

    Returns:
        List[geneset]: The list of genesets.
    """
    
    background_genes = background_geneset.copy()

    with concurrent.futures.ThreadPoolExecutor(max_workers = n_jobs) as executor:

        futures = {executor.submit(_get_multiple_enrichment_data_with_background, 
                                   geneset, 
                                   background_genes, 
                                   database, 
                                   annot_colname, 
                                   annot): geneset for geneset in geneset_list for database in dblist}
        
        pbar = tqdm(total=len(geneset_list), desc='Processing genesets')
        
        for future in as_completed(futures):
            geneset = futures[future]
            pbar.update(1)
        
    pbar.close()

    return geneset_list




# Find terms that contain a given gene

def find_terms_with_gene(gene: str) :

    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/genemap'
    query_string = '?json=true&setup=true&gene=%s'
    gene = gene

    response = requests.get(ENRICHR_URL + query_string % gene)
    
    if not response.ok:
        raise Exception('Error searching for terms')
        
    data = json.loads(response.text)

    return data