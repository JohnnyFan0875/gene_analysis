# Gene_analysis
The script is used for
  1. Gene query: Find genes associated with disease or phenotype from OMIM, ClinVar, Orphanet, and HPO database
  2. Gene name normalization: Normalize gene name from multiple database
  3. Gene depth analysis: Calculate the precentage of gene regions that base coverage < depth


## Requirements
sambamba (gene depth analysis)
pandas
numpy

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
                        genlist file path for analysis. Gene with new line "\n" splitted.
  --output OUTPUT, -o OUTPUT
                        output directory
```

### Output
* <genelist_basename\>_checked.txt: final gene list. filter-out gene will be in the beginning of the file

## Gene depth analysis

```
usage: gene_depth.py [-h] [--depth DEPTH] --output OUTPUT --genelist GENELIST
                     [--bamlist BAMLIST] --dir_bed DIR_BED --mapping_quality
                     {MQ0,MQ10,MQ20,MQ30} [{MQ0,MQ10,MQ20,MQ30} ...]

optional arguments:
  -h, --help            show this help message and exit
  --depth DEPTH, -d DEPTH
                        base depth for analysis. default:20
  --output OUTPUT, -o OUTPUT
                        output directory
  --genelist GENELIST   gene list file. Delimiter "\n"
  --bamlist BAMLIST     bam list file
  --dir_bed DIR_BED, -b DIR_BED
                        directory where bed files locate. bed file name should be <gene>.bed
  --mapping_quality {MQ0,MQ10,MQ20,MQ30} [{MQ0,MQ10,MQ20,MQ30} ...]
```

### Output
* <output\>_<mapping_quality\>.xls




