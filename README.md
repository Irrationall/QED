# QED
Querying  EnrichR Data
- Querying several gene sets to EnrichR at once

<br>

##  Citation
Please cite [EnrichR](https://maayanlab.cloud/Enrichr/)

<br>

## Installation
To install the QED package, follow these steps:

1. Clone the repository (Or download this repository):
    ```
    git clone https://github.com/Irrationall/QED.git
    ```
2. Navigate into the cloned repository:
    ```
    cd QED
    ```
3. Install the package:
    ```
    pip install -e .
    ```

This will install the QED package in editable mode, meaning changes to the source code will be reflected in the installed package without needing to reinstall.

<br>

## Tutorial


```python
from qed.data import readtxt
from qed.data import get_enrichment_dataframes
from qed.data import merge_df
import pandas as pd
```


```python
# Import example data
# Scanpy format: TBD

alist = readtxt("./example_data/endocrinogenesis_rowside.csv", sep=",", format="rowside")
alist = readtxt("./example_data/endocrinogenesis_colside.csv", sep=",", format="colside")

# From 'FindAllMarkers' result from Seurat
alist = readseurat("./example_data/seurat_pbmc_markers.csv", sep=",", index_col=0) # You can pass any pd.DataFrame arguments

```

    Execution time for readtxt: 0.0009999275207519531 seconds
    Execution time for readtxt: 0.0030002593994140625 seconds
    


```python
# Check geneset object
alist
```
* result

    [geneset object [name: Alpha, number of genes: 109, number of GO_dataframes: 0],  
     geneset object [name: Beta, number of genes: 159, number of GO_dataframes: 0],  
     geneset object [name: Delta, number of genes: 51, number of GO_dataframes: 0],  
     geneset object [name: Epsilon, number of genes: 68, number of GO_dataframes: 0],  
     geneset object [name: EP, number of genes: 180, number of GO_dataframes: 0],  
     geneset object [name: Ductal, number of genes: 150, number of GO_dataframes: 0],  
     geneset object [name: Pre-endocrine, number of genes: 61, number of GO_dataframes: 0]]




```python
# Example dblist of EnrichR

dblist = ['KEGG_2021_Human','Reactome_2022','MSigDB_Hallmark_2020','GO_Biological_Process_2021']
```


```python
# Add gene enrichment dataframe to geneset object

aalist = get_enrichment_dataframes(alist, dblist, "Celltype", n_jobs=None)
```


    Processing genesets:   0%|          | 0/7 [00:00<?, ?it/s]



```python
# Let's check geneset object again

aalist
```


* result

    [geneset object [name: Alpha, number of genes: 109, number of GO_dataframes: 4],  
     geneset object [name: Beta, number of genes: 159, number of GO_dataframes: 4],  
     geneset object [name: Delta, number of genes: 51, number of GO_dataframes: 4],  
     geneset object [name: Epsilon, number of genes: 68, number of GO_dataframes: 4],  
     geneset object [name: EP, number of genes: 180, number of GO_dataframes: 4],  
     geneset object [name: Ductal, number of genes: 150, number of GO_dataframes: 4],  
     geneset object [name: Pre-endocrine, number of genes: 61, number of GO_dataframes: 4]]




```python
# mergeing all dataframes at once

df = merge_df(aalist)
df.head(10)
```

| Rank | Term                              | P-value       | Odds ratio | Combined score | Overlapping genes                          | Adjusted p-value | Old p-value | Old adjusted p-value | Database             | Celltype |
|------|-----------------------------------|---------------|------------|----------------|--------------------------------------------|-------------------|-------------|------------------------|----------------------|----------|
| 1    | Pancreas Beta Cells               | 1.888850e-09  | 41.297089  | 829.546906     | [PCSK2, DPP4, SCGN, ABCC8, GCG, IAPP, ISL1] | 4.533240e-08      | 0           | 0                      | MSigDB_Hallmark_2020 | Alpha    |
| 2    | KRAS Signaling Up                 | 1.392889e-05  | 8.126650   | 90.868511      | [RBP4, PCSK1N, TSPAN7, SCG5, USH1C, CPE, SCG3, ... | 1.671466e-04      | 0           | 0                      | MSigDB_Hallmark_2020 | Alpha    |
| 3    | Myogenesis                        | 2.394210e-02  | 3.827988   | 14.286499      | [CAMK2B, NQO1, GPX3, NCAM1]                | 1.758437e-01      | 0           | 0                      | MSigDB_Hallmark_2020 | Alpha    |
| 4    | Reactive Oxygen Species Pathway   | 2.930728e-02  | 7.891827   | 27.857515      | [NQO1, GPX3]                               | 1.758437e-01      | 0           | 0                      | MSigDB_Hallmark_2020 | Alpha    |
| 5    | Spermatogenesis                   | 3.763212e-02  | 4.236492   | 13.895260      | [PCSK1N, SCG5, SCG3]                       | 1.806342e-01      | 0           | 0                      | MSigDB_Hallmark_2020 | Alpha    |
| 6    | Epithelial Mesenchymal Transition | 9.615827e-02  | 2.829327   | 6.625603       | [RGS4, CDH2, SCG2]                         | 3.310276e-01      | 0           | 0                      | MSigDB_Hallmark_2020 | Alpha    |
| 7    | Protein Secretion                  | 9.655232e-02  | 3.936568   | 9.202398       | [OCRL, PAM]                               | 3.310276e-01      | 0           | 0                      | MSigDB_Hallmark_2020 | Alpha    |
| 8    | Peroxisome                         | 1.103425e-01  | 3.626351   | 7.993080       | [TTR, ABCC8]                              | 3.310276e-01      | 0           | 0                      | MSigDB_Hallmark_2020 | Alpha    |
| 9    | Hedgehog Signaling                 | 1.787346e-01  | 5.252910   | 9.044741       | [SCG2]                                    | 4.766255e-01      | 0           | 0                      | MSigDB_Hallmark_2020 | Alpha    |
| 10   | Xenobiotic Metabolism              | 2.975500e-01  | 1.859058   | 2.253500       | [NQO1, RBP4]                              | 6.491999e-01      | 0           | 0                      | MSigDB_Hallmark_2020 | Alpha    |

