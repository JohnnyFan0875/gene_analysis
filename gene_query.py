###find gene by keyword
import os,sys
import time
import xml.etree.ElementTree as ET
import argparse

def get_parser():

    parser = argparse.ArgumentParser()
    parser.add_argument('--keyword','-k',required=True,type=str,help='keyword to query in database. e.g. "Kennedy disease"')
    parser.add_argument('--output','-o',required=True,type=str,help='output directory')
    args = parser.parse_args()

    args.output = os.path.abspath(args.output.strip())
    keyword_li = args.keyword.lower().split(' ')

    return args, keyword_li

class db_download():

    def hpo_db(self):

        if 'g2p_'+time_Ym+'.xls' in os.listdir(os.path.join(script_path,'database')):
            print('Latest HPO database already in folder. Skip')
        else: 
            os.system('wget -c -q --show-progress http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt -O ' + os.path.join(script_path,'database',f'g2p_{time_Ym}.xls'))
            print('HPO database downloaded.')

    def clinvar_db(self):

        if 'clinvar_'+time_Ym+'.txt' in os.listdir(os.path.join(script_path,'database')):
            print('Latest ClinVar database already in folder. Skip')
        else:
            os.system('wget -c -q --show-progress https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz -O ' + os.path.join(script_path,'database',f'clinvar_{time_Ym}.txt.gz'))
            os.system('gunzip ' + os.path.join(script_path,'database',f'clinvar_{time_Ym}.txt.gz'))
            print('ClinVar database downloaded.')

    def orphanet_db(self):

        if 'orphanet_'+time_Ym+'.xls' in os.listdir(os.path.join(script_path,'database')):
            print('Latest Orphanet database already in folder. Skip')
        else:
            os.system('wget -c -q --show-progress https://www.orphadata.com/data/xml/en_product6.xml -O ' + os.path.join(script_path,'database',f'orphanet_{time_Ym}.xml'))
            
            with open(os.path.join(script_path,'database','orphanet_'+time_Ym+'.xls'),'w',encoding='utf-8') as f_output:

                f_output.write('Orpha_code\tOrpha_name\tGene\n')
            
                tree = ET.parse(os.path.join(script_path,'database','orphanet_{}.xml'.format(time_Ym)))
                root = tree.getroot()

                for disorder in root.find('DisorderList').findall('Disorder'):
                    code = disorder.find('OrphaCode').text
                    name = disorder.find('Name').text
                    li_gene = list()
                    gene_tags = disorder.find('DisorderGeneAssociationList').findall('DisorderGeneAssociation')
                    for gene_tag in gene_tags:
                        li_gene.append(gene_tag.find('Gene').find('Symbol').text)

                    f_output.write('\t'.join([code,name,','.join(li_gene)]) + '\n')

class db_info_query():

    def hpo_g2p(self):

        g2p_fp = os.path.join(script_path,'database',sorted([i for i in os.listdir(os.path.join(script_path,'database')) if i.startswith('g2p')])[-1])
        hpo_g2p_li = list()

        with open(g2p_fp,'r',encoding='utf-8') as g2p_input:
            g2p_line = g2p_input.readline()
            while g2p_line:
                if not '#' in g2p_line:
                    gene, hpo_id, hpo_name = [g2p_line.split('\t')[i].strip() for i in range(1,4)]
                    for disorder in [hpo_name.lower()]:
                        if all(key in disorder for key in keyword_li):
                            hpo_g2p_li.append(['HPO','.','.',hpo_id,hpo_name,gene])
                            break
                g2p_line = g2p_input.readline()

        print('HPO g2p file query completed.')
        return hpo_g2p_li

    def clinvar(self):

        clinvar_fp = os.path.join(script_path,'database',sorted([i for i in os.listdir(os.path.join(script_path,'database')) if i.startswith('clinvar')])[-1])
        clinvar_li = list()
        with open(clinvar_fp,'r',encoding='utf-8') as clinvar_input:
            clinvar_line = clinvar_input.readline()
            
            while clinvar_line:

                if clinvar_line.startswith('#'):
                    gene_idx = clinvar_line.split('\t').index('Name')
                    variant_id_idx = clinvar_line.split('\t').index('VariationID')
                    phenotype_idx = clinvar_line.split('\t').index('PhenotypeList')
                    pathogenicity_idx = clinvar_line.split('\t').index('ClinicalSignificance')
                    pass

                else:
                    try:
                        clinvar_variant = clinvar_line.split('\t')[gene_idx]
                        clinvar_gene = clinvar_line.split('\t')[gene_idx].split('(')[1].split(')')[0].strip()
                        clinvar_ID = clinvar_line.split('\t')[variant_id_idx].strip()
                        clinvar_phenotype = clinvar_line.split('\t')[phenotype_idx].strip().lower()
                        clinvar_pathogenicity = clinvar_line.split('\t')[pathogenicity_idx].strip().lower()

                        if ('?' in clinvar_gene) or ('chr' in clinvar_gene and ':' in clinvar_gene and '-' in clinvar_gene):
                            pass
                        elif all(i in clinvar_variant for i in ['g.(',')_(']):
                            pass
                        elif all(i in clinvar_variant for i in ['c.(',')_(']):
                            pass
                        elif 'pathogenic' not in clinvar_pathogenicity or 'conflict' in clinvar_pathogenicity:
                            pass
                        else:
                            for disorder in clinvar_phenotype.split('|'):
                                if all(key in disorder for key in keyword_li):
                                    clinvar_li.append(['ClinVar',clinvar_ID,clinvar_phenotype,'.','.',clinvar_gene])
                                    break
                    except:
                        pass 

                clinvar_line = clinvar_input.readline()

        print('ClinVar file query completed.')
        return clinvar_li

    def omim(self):

        omim_fp = os.path.join(script_path,'database',sorted([i for i in os.listdir(os.path.join(script_path,'database')) if i.startswith('omim')])[-1])
        omim_li = list()

        with open(omim_fp,'r',encoding='utf-8') as omim_input:
            omim_line = omim_input.readline()
            while omim_line:
                if not '#' in omim_line:
                    omim_gene = omim_line.split('\t')[0].strip()
                    omim_name = omim_line.split('\t')[1].strip().lower()

                    for disorder in [omim_name]:
                        if all(key in disorder for key in keyword_li):
                            omim_li.append(['OMIM','.',omim_name,'.','.',omim_gene])
                            break
                omim_line = omim_input.readline()

        print('OMIM file query completed.')
        return omim_li

    def orphanet(self):

        orpha_fp = os.path.join(script_path,'database',sorted([i for i in os.listdir(os.path.join(script_path,'database')) if i.startswith('orphanet') and i.endswith('.xls')])[-1])
        orphanet_li = list()

        with open(orpha_fp,'r',encoding='utf-8') as orpha_input:
            orpha_line = orpha_input.readline()
            while orpha_line:
                if not '#' in orpha_line:
                    orpha_genes = orpha_line.split('\t')[2].strip()
                    orpha_name = orpha_line.split('\t')[1].lower().strip()

                    for disorder in [orpha_name]:
                        if all(key in disorder for key in keyword_li):
                            for gene in orpha_genes.split(','):
                                orphanet_li.append(['Orphanet','.',orpha_name,'.','.',gene])
                            break
                orpha_line = orpha_input.readline()

        print('Orphanet file query completed.')
        return orphanet_li

def merge_db_query(query_db):

    os.makedirs(args.output,exist_ok=True)

    genelist = list()
    output_detail_fp = os.path.join(args.output,'detail_' + ''.join(keyword_li) + '.xls')
    output_genelist_fp = os.path.join(args.output,'genelist_' + ''.join(keyword_li) + '.txt')

    with open(output_detail_fp,'w',encoding='utf-8') as output_detail_file:

        output_detail_file.write(f'Source\tDisease_ID\tDisease_Name\tHPO_ID\tHPO_Name\tGene\n')

        d = {'HPO':query_db[0],'OMIM':query_db[1],'Orphanet':query_db[2],'ClinVar':query_db[3]}
        
        for db,db_li in d.items():
            db_genelist = sorted(dict.fromkeys([i[5] for i in db_li]))
            output_detail_file.write(f'{db}\tTotal\t\t\t\t' + ','.join(db_genelist) + '\n')

            for line in db_li:
                output_detail_file.write('\t'.join(line)+'\n')
                genelist += line[5].strip().split(',')

    total_genelist = sorted(dict.fromkeys(genelist))
    with open(output_genelist_fp,'w',encoding='utf-8') as output_genelist_file:
        output_genelist_file.write('\n'.join(total_genelist))

    print('detail and genelist files generated.')



if __name__ == "__main__":
    
    start_time = time.time()

    script_path = os.path.dirname(os.path.abspath(__file__))
    time_Ym = time.strftime("%Y%m")

    args, keyword_li = get_parser()

    #download database
    print('Database downloading...')
    db_download = db_download()
    db_download.hpo_db()
    db_download.clinvar_db()
    db_download.orphanet_db()

    #database query
    print('Gene query from database...')
    db_info_query = db_info_query()
    hpo_g2p_li = db_info_query.hpo_g2p()
    clinvar_li = db_info_query.clinvar()
    omim_li = db_info_query.omim()
    orphanet_li = db_info_query.orphanet()
    
    #make final result file
    print('Generate final result...')
    query_db = (hpo_g2p_li,omim_li,orphanet_li,clinvar_li)
    query_genelist = merge_db_query(query_db)

    end_time = time.time()
    print('*'*30+'\nProgram finished. Total spent time %0.2fs for this program.'%(end_time-start_time))