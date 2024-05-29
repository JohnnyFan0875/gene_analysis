# Gene_analysis
The script is used for
  1. Gene query: Find genes associated with disease or phenotype from OMIM, ClinVar, Orphanet, and HPO database.
  2. Gene name normalization: Normalize gene name from multiple database.
  3. Gene depth analysis: Calculate the coverage of the regions of genes. 

## Gene query
```
$ python3 gene_query.py --keyword keyword --output output_folder_name
```
Arguments
* --keyword: Keyword of disease or phenotype. Double quotes ("") is required if keyword contains space.
* --output: Output folder name. The folder will be created in ./data

Note
* ClinVar, Orphanet and HPO database will be downloaded automatically.
* Latest OMIM data should be acquired from the website (https://www.omim.org/downloads).
  OMIM database (2021/07) is in ./database
