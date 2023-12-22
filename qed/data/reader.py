from .structure import geneset
from typing import List, Optional
import pandas as pd
import time


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
            genes = list(filter(None, line.split(",")[1:]))
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
def readseurat() -> List[geneset] :
    
    setlist = []
    
    return setlist


@calc_time
def readscanpy() -> List[geneset] :
    
    setlist = []
    
    return setlist