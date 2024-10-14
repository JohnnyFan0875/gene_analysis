import os
import sys
import argparse


    
def get_parser():

    parser = argparse.ArgumentParser()
    parser.add_argument('--variant_level',type=int,help='Review status of variant in ClinVar. Genes are collected if there is at least one variant present at the specific variant level or above.')
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--keyword_filepath',type=str,help='Path to the text file containing keywords (one per line) for finding associated genes.')
    group.add_argument('--genelist_filepath',type=str,help='Path to the text file containing genes (one per line) for selection based on variant level and print out variant phenotype.')

    args = parser.parse_args()
    return args



def keyword_associated_gene():

    keyword_filepath = os.path.abspath(args.keyword_filepath.strip())
    with open(keyword_filepath,'r',encoding='utf-8') as f_input:
        keyword_d = {i.strip():{'all':set(),'certificate':[]} for i in f_input.read().strip().split('\n')}

    #clinvar
    with open(clinvar_filepath,'r', encoding='utf-8') as f_input:
        line = f_input.readline()
        while line:

            gene_check = False
            pathogenicity = line.split('\t')[5].lower()
            annotation = line.split('\t')[6]
            phenotype = line.split('\t')[9].lower()
            status = line.split('\t')[10].strip().replace('☆☆☆☆','0').replace('★☆☆☆','1').replace('★★☆☆','2').replace('★★★☆','3').replace('★★★★','4')

            if '(' in annotation and ')' in annotation and '?' not in annotation:
                gene = annotation.split('(')[1].split(')')[0]
                gene_check = True
            else:
                pass

            if gene_check and 'pathogenic' in pathogenicity and 'conflict' not in pathogenicity and int(status)>=args.variant_level:
                for keyword in keyword_d:
                    for single_pheno in phenotype.split('|'):
                        if all(i.lower() in single_pheno.lower() for i in keyword.split(' ')):
                            keyword_d[keyword]['all'].add(gene)
                            break
                            
            line = f_input.readline()

    #output
    output_filepath = keyword_filepath + '_final.xls'
    with open(output_filepath,'w',encoding='utf-8') as f_output:
        f_output.write('\t'.join(['Keyword','All_gene_no','All_gene_list','Certificate_gene_no','Certificate_gene_list']) + '\n')
        for keyword in keyword_d:
            for gene in keyword_d[keyword]['all']:
                if gene in certificate_li:
                    keyword_d[keyword]['certificate'].append(gene)
                else:
                    pass
            f_output.write('\t'.join([keyword,
                                      str(len(list(keyword_d[keyword]['all']))),
                                      ','.join(sorted(list(keyword_d[keyword]['all']))),
                                      str(len(keyword_d[keyword]['certificate'])),
                                      ','.join(sorted(list(keyword_d[keyword]['certificate'])))
                                      ]) + '\n')
            


def get_gene_clinvar_phenotype():

    genelist_filepath = os.path.abspath(args.genelist_filepath.strip())
    with open(genelist_filepath,'r',encoding='utf-8') as f_input:
        gene_dict = {i:{} for i in f_input.read().strip().split('\n')}

    #get phenotype in clinvar
    with open(clinvar_filepath,'r', encoding='utf-8') as f_input:
        line = f_input.readline()
        while line:

            gene_check = False
            pathogenicity = line.split('\t')[5].lower()
            annotation = line.split('\t')[6]
            phenotype = line.split('\t')[9].lower()
            status = line.split('\t')[10].strip().replace('☆☆☆☆','0').replace('★☆☆☆','1').replace('★★☆☆','2').replace('★★★☆','3').replace('★★★★','4')

            if '(' in annotation and ')' in annotation and '?' not in annotation:
                gene = annotation.split('(')[1].split(')')[0]
                gene_check = True
            else:
                pass

            if gene_check and gene in gene_dict and 'pathogenic' in pathogenicity and 'conflict' not in pathogenicity and int(status)>=args.variant_level:
                gene_dict[gene][annotation] = phenotype
                            
            line = f_input.readline() 

    for gene in gene_dict:
        if gene_dict[gene] == {}:
            print(f'{gene} has no {args.variant_level}+star pathogenic variant in clinvar')         

    #output
    output_filepath = genelist_filepath + '_final.xls'
    with open(output_filepath,'w',encoding='utf-8') as f_output:
        f_output.write('\t'.join(['Gene','Annotation','Phenotype']) + '\n')
        for gene in gene_dict:
            for annotation in gene_dict[gene]:
                f_output.write('\t'.join([gene,annotation,gene_dict[gene][annotation]]) + '\n')

        

if __name__ == "__main__":
    
    script_path = os.path.dirname(os.path.abspath(__file__))
    
    args = get_parser()

    clinvar_filepath = sorted([os.path.join('/home/khhg/LVA/Database/',i) for i in os.listdir('/home/khhg/LVA/Database/') if i.startswith('clinvar_') and i.endswith('.txt')])[-1]

    with open(os.path.join(script_path,'WES_certificate_gene.txt'),'r',encoding='utf-8') as f_input:
        certificate_li = f_input.read().strip().split(',')

    if vars(args)['keyword_filepath'] is not None:
        keyword_associated_gene()
    else:
        get_gene_clinvar_phenotype()