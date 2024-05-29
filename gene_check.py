import os,sys
import time

class parser():

    def parser(self):

        import argparse

        descrip = 'Gene check from database. Command line: python3 /home/johnny/gene_query_check_depth/gene_check_v2.py'
        parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description=descrip)
        parser.add_argument('--genelist', type=str, required=True, help='Required. Genlist file path for analysis. Gene with new line "\\n" splitted.')
        args = parser.parse_args()
        return args

class gene_check():

    #get refGene and chrM_ensembl information
    def db_info_compile(self):

        All_genelist = set()
        with open(os.path.join(script_path,'database','refGene_20220505.txt'), 'r') as refgene_input,\
             open(os.path.join(script_path,'database','chrM_ensembl.tsv'),'r') as chrM_ensembl_input:

            for row in refgene_input.readlines():
                if row.split('\t')[1].startswith('NM_') and '_' not in row.split('\t')[2]:
                    All_genelist.add(row.split('\t')[12].strip())

            for row in chrM_ensembl_input.readlines():
                All_genelist.add(row.strip().split('\t')[-1].strip())
        
        return All_genelist

    #replace new gene name with old gene name
    def gene_name_update(self,genelist):
        
        with open(os.path.join(script_path,'database','update_manual.txt'),'r') as update_manual_input,\
             open(os.path.join(script_path,'database','chrM_alias.txt'),'r') as chrM_alias_input:

            for chrM_alias in chrM_alias_input.readlines():
                if chrM_alias.split('\t')[1].strip() in genelist:
                    genelist.remove(chrM_alias.split('\t')[1].strip())
                    genelist.append(chrM_alias.split('\t')[0].strip())
                    print(chrM_alias.split('\t')[1].strip() + '-->' + chrM_alias.split('\t')[0].strip())


            for update_row in update_manual_input.readlines():
                if update_row.split('\t')[0].strip() in genelist:
                    genelist.remove(update_row.split('\t')[0].strip())
                    genelist.append(update_row.split('\t')[1].strip())
                    print(update_row.split('\t')[0].strip() + '-->' + update_row.split('\t')[1].strip())

        print('Gene name renewed.')
        return genelist

    #remove duplicate gene
    def gene_rmdup(self,genelist):

        genelist2 = sorted(list(dict.fromkeys(genelist)))
        if len(genelist2) != len(genelist):
            print(f'Gene list number from {len(genelist)} to {len(genelist2)} via rmdup step.')
        else:
            print('No duplicate gene found.')

        return genelist2

    #gene classification
    def gene_classification(self, genelist, All_genelist):

        #get genecard database information
        genecard_d = dict()
        with open(os.path.join(script_path,'database','genecard_20220725.txt'),'r') as genecard_input:
            for row in genecard_input.readlines()[1:]:
                genecard_d[row.split('\t')[0].strip()] = row.split('\t')[1].strip()

        #define category
        classification_d = {'final':[list(),'Final_genelist'],
        'filout1':[list(),'### Not in refgene db, genecard defined this as protein coding gene'],
        'filout2':[list(),'### Not in refgene db, mitochondria gene'],
        'filout3':[list(),'### Not in refgene db, not protein coding gene or MT gene'],
        'filout4':[list(),'### Not in refgene db and genecard db']}
    
        for gene in genelist:
            if gene in All_genelist:
                classification_d['final'][0] += [gene]
            else:
                if gene in genecard_d:
                    if genecard_d[gene] == 'Protein Coding':
                        classification_d['filout1'][0] += [gene]
                    elif gene.startswith('MT-'):
                        classification_d['filout2'][0] += [gene]
                    else:
                        classification_d['filout3'][0] += [gene]
                else:
                    classification_d['filout4'][0] += [gene]

        return classification_d

    #output gene check information file
    def output_gene_check(self, classification_d, file_path):

        output_fp = file_path.replace('.xls','').replace('.txt','')+'_checked.txt'
        with open(output_fp,'w',encoding='utf-8') as output_file:
            for filout in [i for i in classification_d if i != 'final']:
                for gene in classification_d[filout][0]:
                    output_file.write(gene+'\t'+classification_d[filout][1]+'\n')
                    print(gene+'\t'+classification_d[filout][1])
            output_file.write('\n'.join(classification_d['final'][0]))

        return classification_d['final'][0]

if __name__ == "__main__":

    start_time = time.time()
    parser = parser()
    script_path = os.path.dirname(os.path.abspath(__file__))
    
    file_path = parser.parser().genelist
    if not os.path.isfile(os.path.abspath(file_path.strip())):
        sys.exit('Requested sample file path does not exist. Program stop.')
    with open(file_path.strip(),'r',encoding='utf-8') as genelist_file_input:
        genelist = [i.strip() for i in genelist_file_input.read().strip().split('\n')]  

    gene_check = gene_check()
    All_genelist = gene_check.db_info_compile()
    genelist = gene_check.gene_name_update(genelist)
    genelist = gene_check.gene_rmdup(genelist)
    classification_d = gene_check.gene_classification(genelist, All_genelist)
    checked_genelist = gene_check.output_gene_check(classification_d, file_path)

    end_time = time.time()
    print('*'*30+'\nProgram finished. Total spent time {}s for this program.'.format(end_time-start_time))
    