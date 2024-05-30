import os,sys
import time
import argparse

def get_parser():

    parser = argparse.ArgumentParser()
    parser.add_argument('--genelist', '-g', type=str, required=True, help='genlist file path for analysis. Gene with new line "\\n" splitted.')
    parser.add_argument('--output','-o',required=True,type=str,help='output directory')
    args = parser.parse_args()

    args.genelist = os.path.abspath(args.genelist.strip())
    args.output = os.path.abspath(args.output.strip())

    if not os.path.isfile(args.genelist):
        sys.exit('Gene list file not found. Program stop.')

    return args

def get_genelist():

    with open(args.genelist,'r',encoding='utf-8') as genelist_file_input:
        genelist = [i.strip() for i in genelist_file_input.read().strip().split('\n')]  

    return genelist

class gene_check():

    #get refGene and chrM_ensembl information
    def db_info_compile(self):

        print('Get refGene and ensembl_chrM information...')

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

        print('Update gene name...')
        
        with open(os.path.join(script_path,'database','update_manual.txt'),'r') as update_manual_input,\
             open(os.path.join(script_path,'database','chrM_alias.txt'),'r') as chrM_alias_input:

            for chrM_alias in chrM_alias_input.readlines():
                if chrM_alias.split('\t')[1].strip() in genelist:
                    genelist.remove(chrM_alias.split('\t')[1].strip())
                    genelist.append(chrM_alias.split('\t')[0].strip())
                    print('*'*10 + chrM_alias.split('\t')[1].strip() + '-->' + chrM_alias.split('\t')[0].strip())

            for update_row in update_manual_input.readlines():
                if update_row.split('\t')[0].strip() in genelist:
                    genelist.remove(update_row.split('\t')[0].strip())
                    genelist.append(update_row.split('\t')[1].strip())
                    print('*'*10 + update_row.split('\t')[0].strip() + '-->' + update_row.split('\t')[1].strip())

        return genelist

    #remove duplicate gene
    def gene_rmdup(self,genelist):

        print('Gene duplicate removal...')
        
        genelist2 = sorted(list(dict.fromkeys(genelist)))
        if len(genelist2) != len(genelist):
            print(f'Gene list number from {len(genelist)} to {len(genelist2)} via rmdup step.')
        else:
            print('No duplicate gene found.')

        return genelist2

    #gene classification
    def gene_classification(self,genelist):

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
    def output_gene_check(self,genelist):

        os.makedirs(args.output,exist_ok=True)
        output_fp = os.path.join(
            args.output,
            os.path.splitext(os.path.basename(args.genelist))[0] + '_checked.txt'
        )

        print('Final gene number:' + str(len(classification_d['final'][0])))
        print('Final gene list: ' + ','.join(classification_d['final'][0]))

        with open(output_fp,'w',encoding='utf-8') as output_file:
            for filout in [i for i in classification_d if i != 'final']:
                for gene in classification_d[filout][0]:
                    output_file.write(gene+'\t'+classification_d[filout][1]+'\n')
                    print(gene+'\t'+classification_d[filout][1])
            output_file.write('\n'.join(classification_d['final'][0]))



if __name__ == "__main__":

    start_time = time.time()

    args = get_parser()
    script_path = os.path.dirname(os.path.abspath(__file__))

    #get gene list
    genelist = get_genelist()
    
    #gene check
    gene_check = gene_check()
    All_genelist = gene_check.db_info_compile()
    genelist = gene_check.gene_name_update(genelist)
    genelist = gene_check.gene_rmdup(genelist)

    #gene classification
    classification_d = gene_check.gene_classification(genelist)
    checked_genelist = gene_check.output_gene_check(genelist)

    end_time = time.time()
    print('*'*30+'\nProgram finished. Total spent time %0.2fs for this program.'%(end_time-start_time))
    