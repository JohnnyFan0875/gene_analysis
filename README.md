# Gene_analysis
The script is used for
  1. Gene query: Find genes associated with disease or phenotype from OMIM, ClinVar, Orphanet, and HPO database.
  2. Gene name normalization: Normalize gene name from multiple database.
  3. Gene depth analysis: Calculate the coverage of the regions of genes. 

## Gene query
```
usage: gene_query.py [-h] --keyword KEYWORD --output OUTPUT

Gene query from database. Command line: python3 gene_query.py

optional arguments:
  -h, --help            show this help message and exit
  --keyword KEYWORD, -k KEYWORD
                        keyword to query in database. e.g. "Kennedy disease"
  --output OUTPUT, -o OUTPUT
                        output directory
```

Note
* Latest OMIM data should be acquired from the website (https://www.omim.org/downloads).  
  OMIM database (2021/07) already in database

### Output
* detail_<keyword>.xls: gene query detail in every database
* genelist_<keyword>.txt: total gene from all database 

