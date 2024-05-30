# Gene_analysis
The script is used for
  1. Gene query: Find genes associated with disease or phenotype from OMIM, ClinVar, Orphanet, and HPO database
  2. Gene name normalization: Normalize gene name from multiple database
  3. Gene depth analysis: Calculate the coverage of the regions of genes

## Gene query
```
usage: gene_query.py [-h] --keyword KEYWORD --output OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  --keyword KEYWORD, -k KEYWORD
                        keyword to query in database. e.g. "Kennedy disease"
  --output OUTPUT, -o OUTPUT
                        output directory
```

Note
* Latest OMIM data should be acquired from the website (https://www.omim.org/downloads) 
  OMIM database (2021/07) in database

### Output
* detail_<keyword\>.xls: detail gene information from all database
* genelist_<keyword\>.txt: total gene from all database 

## Gene name normalization
```
usage: gene_check.py [-h] --genelist GENELIST --output OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  --genelist GENELIST, -g GENELIST
                        genlist file path for analysis. Gene with new line
                        "\n" splitted.
  --output OUTPUT, -o OUTPUT
                        output directory
```

### Output
* <genelist_basename\>_checked.txt: final gene list. filter-out gene will be in the beginning of the file

