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
from qed.data.reader import readtxt
from qed.data.query import get_enrichment_dataframes
from qed.data.structure import merge_df
import pandas as pd
```


```python
# Import example data
# Seurat and Scanpt format: TBD

alist = readtxt("./example_data/endocrinogenesis_rowside.csv", sep=",", format="rowside")
blist = readtxt("./example_data/endocrinogenesis_colside.csv", sep=",", format="colside")
```

    Execution time for readtxt: 0.0009999275207519531 seconds
    Execution time for readtxt: 0.0030002593994140625 seconds
    


```python
# Check geneset object
alist
```




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



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Rank</th>
      <th>Term</th>
      <th>P-value</th>
      <th>Odds ratio</th>
      <th>Combined score</th>
      <th>Overlapping genes</th>
      <th>Adjusted p-value</th>
      <th>Old p-value</th>
      <th>Old adjusted p-value</th>
      <th>Database</th>
      <th>Celltype</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>Pancreas Beta Cells</td>
      <td>1.888850e-09</td>
      <td>41.297089</td>
      <td>829.546906</td>
      <td>[PCSK2, DPP4, SCGN, ABCC8, GCG, IAPP, ISL1]</td>
      <td>4.533240e-08</td>
      <td>0</td>
      <td>0</td>
      <td>MSigDB_Hallmark_2020</td>
      <td>Alpha</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2</td>
      <td>KRAS Signaling Up</td>
      <td>1.392889e-05</td>
      <td>8.126650</td>
      <td>90.868511</td>
      <td>[RBP4, PCSK1N, TSPAN7, SCG5, USH1C, CPE, SCG3,...</td>
      <td>1.671466e-04</td>
      <td>0</td>
      <td>0</td>
      <td>MSigDB_Hallmark_2020</td>
      <td>Alpha</td>
    </tr>
    <tr>
      <th>2</th>
      <td>3</td>
      <td>Myogenesis</td>
      <td>2.394210e-02</td>
      <td>3.827988</td>
      <td>14.286499</td>
      <td>[CAMK2B, NQO1, GPX3, NCAM1]</td>
      <td>1.758437e-01</td>
      <td>0</td>
      <td>0</td>
      <td>MSigDB_Hallmark_2020</td>
      <td>Alpha</td>
    </tr>
    <tr>
      <th>3</th>
      <td>4</td>
      <td>Reactive Oxygen Species Pathway</td>
      <td>2.930728e-02</td>
      <td>7.891827</td>
      <td>27.857515</td>
      <td>[NQO1, GPX3]</td>
      <td>1.758437e-01</td>
      <td>0</td>
      <td>0</td>
      <td>MSigDB_Hallmark_2020</td>
      <td>Alpha</td>
    </tr>
    <tr>
      <th>4</th>
      <td>5</td>
      <td>Spermatogenesis</td>
      <td>3.763212e-02</td>
      <td>4.236492</td>
      <td>13.895260</td>
      <td>[PCSK1N, SCG5, SCG3]</td>
      <td>1.806342e-01</td>
      <td>0</td>
      <td>0</td>
      <td>MSigDB_Hallmark_2020</td>
      <td>Alpha</td>
    </tr>
    <tr>
      <th>5</th>
      <td>6</td>
      <td>Epithelial Mesenchymal Transition</td>
      <td>9.615827e-02</td>
      <td>2.829327</td>
      <td>6.625603</td>
      <td>[RGS4, CDH2, SCG2]</td>
      <td>3.310276e-01</td>
      <td>0</td>
      <td>0</td>
      <td>MSigDB_Hallmark_2020</td>
      <td>Alpha</td>
    </tr>
    <tr>
      <th>6</th>
      <td>7</td>
      <td>Protein Secretion</td>
      <td>9.655232e-02</td>
      <td>3.936568</td>
      <td>9.202398</td>
      <td>[OCRL, PAM]</td>
      <td>3.310276e-01</td>
      <td>0</td>
      <td>0</td>
      <td>MSigDB_Hallmark_2020</td>
      <td>Alpha</td>
    </tr>
    <tr>
      <th>7</th>
      <td>8</td>
      <td>Pperoxisome</td>
      <td>1.103425e-01</td>
      <td>3.626351</td>
      <td>7.993080</td>
      <td>[TTR, ABCC8]</td>
      <td>3.310276e-01</td>
      <td>0</td>
      <td>0</td>
      <td>MSigDB_Hallmark_2020</td>
      <td>Alpha</td>
    </tr>
    <tr>
      <th>8</th>
      <td>9</td>
      <td>Hedgehog Signaling</td>
      <td>1.787346e-01</td>
      <td>5.252910</td>
      <td>9.044741</td>
      <td>[SCG2]</td>
      <td>4.766255e-01</td>
      <td>0</td>
      <td>0</td>
      <td>MSigDB_Hallmark_2020</td>
      <td>Alpha</td>
    </tr>
    <tr>
      <th>9</th>
      <td>10</td>
      <td>Xenobiotic Metabolism</td>
      <td>2.975500e-01</td>
      <td>1.859058</td>
      <td>2.253500</td>
      <td>[NQO1, RBP4]</td>
      <td>6.491999e-01</td>
      <td>0</td>
      <td>0</td>
      <td>MSigDB_Hallmark_2020</td>
      <td>Alpha</td>
    </tr>
  </tbody>
</table>
</div>


