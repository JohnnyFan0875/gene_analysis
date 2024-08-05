import os
import sys
import argparse
import gzip

def get_parser():

    parser = argparse.ArgumentParser()
    parser.add_argument('--genelist','-g',type=str,required=True,help='file with gene splitted by "\\n"')
    args = parser.parse_args()

    return args

def get_genelist():

    args.genelist = os.path.abspath(args.genelist)
    if not os.path.exists(args.genelist):
        sys.exit('genelist file not found. Program stop.')
    
    with open(args.genelist,'r',encoding='utf-8') as genelist_input:
        genelist = genelist_input.read().strip().split('\n')

    for gene in genelist: #check if bed files already exist
        if os.path.exists(os.path.join(script_path,'database','gene_bed','MGIcapture',f'{gene}.bed')) and\
           os.path.exists(os.path.join(script_path,'database','gene_bed','cdsplice5',f'{gene}.bed')):
           genelist.remove(gene)

    return genelist

class create_bed():

    def check_gene(func):
        def wrap(self,genelist):
            bed_d, bed_abbr = func(self,genelist)
            check_li = list()
            if bed_abbr == 'Mitochondrial':
                for gene in [i for i in genelist if i.startswith('MT-')]:
                    if bed_d[gene] == []:
                        check_li.append(gene)
            else:
                for gene in [i for i in genelist if not i.startswith('MT-')]:
                    if bed_d[gene] == []:
                        check_li.append(gene)
            if check_li != []:
                sys.exit(', '.join(check_li) + ' gene(s) not found in ' + bed_abbr + ' bed file. Program stop.')
            print(f'{bed_abbr} bed file dictionary gene check complete.')
            return bed_d
        return wrap

    def prepare_dict(self):

        update_d = {} #update_d[new_gene_name] = old_gene_name
        transcript_d = {} #transcript_d[gene] = transcript
        ncbi_d = {} #ncbi_d[new] = [old_gene_name]

        with open(os.path.join(script_path,'database','update_manual.txt'),'r',encoding='utf-8') as update_input:
            for line in update_input.read().strip().split('\n'):
                update_d[line.split('\t')[1]] = line.split('\t')[0]

        with open(os.path.join(script_path,'database','canonical_refseq.tsv'),'r',encoding='utf-8') as transcript_input:
            for line in transcript_input.readlines():
                transcript_d[line.split('\t')[3]] =line.split('\t')[4].split('.')[0]

        with gzip.open(os.path.join(script_path,'database','Homo_sapiens.gene_info.gz'),'r') as ncbi_input:
            line = ncbi_input.readline()
            while line:
                line = line.decode()
                if not line.split('\t')[4] == '-':
                    ncbi_d[line.split('\t')[2].strip()] = line.split('\t')[4].strip().split('|')
                line = ncbi_input.readline()

        return update_d, transcript_d, ncbi_d

    @check_gene
    def cdsplice5(self,genelist):

        bed_path = '/home/khhg/LVA/Database/MGI_cds_splice5.bed'
        bed_d = {}
        for gene in [i for i in genelist if not i.startswith('MT-')]:
            bed_d[gene] = list()
            with open(bed_path,'r',encoding='utf-8') as bed_input:
                line = bed_input.readline()
                while line:
                    input_gene = line.split('\t')[3].strip()
                    if input_gene == gene:
                        bed_d[gene] += [line.strip()]
                    elif gene in update_d and input_gene == update_d[gene]: #bed file may use old gene name
                        bed_d[gene] += [line.strip()]
                    line = bed_input.readline()

        print('cdsplice5 bed file dictionary created.')
        return bed_d, 'cdsplice5'

    @check_gene
    def MGIcapture(self,genelist):

        bed_path = '/home/khhg/LVA/Database/MGI_Exome_Capture_V5.anotation.bed'
        bed_d = {}
        for gene in [i for i in genelist if not i.startswith('MT-')]:
            bed_d[gene] = list()
            with open(bed_path,'r',encoding='utf-8') as bed_input:
                line = bed_input.readline()
                while line:
                    content_li = line.strip().split('\t')[3].replace('|',',').split(',')
                    if gene in content_li:
                        bed_d[gene] += [line.strip()] 
                    elif gene in update_d and update_d[gene] in content_li: #bed file may use old gene name
                        bed_d[gene] += [line.strip()]
                    elif gene in transcript_d and transcript_d[gene] in content_li:
                        bed_d[gene] += [line.strip()]
                    elif gene in update_d and update_d[gene] in transcript_d and transcript_d[update_d[gene]] in content_li:
                        bed_d[gene] += [line.strip()]
                    elif gene in ncbi_d and any(i in content_li for i in ncbi_d[gene]):
                        bed_d[gene] += [line.strip()]
                    line = bed_input.readline()

        print('MGIcapture bed file dictionary created.')
        return bed_d, 'MGIcapture'

    @check_gene
    def mitochondria(self,genelist):

        bed_path = '/home/khhg/LVA/Database/chrM_ensembl.tsv'
        bed_d = {}
        for gene in [i for i in genelist if i.startswith('MT-')]:
            bed_d[gene] = list()
            with open(bed_path,'r',encoding='utf-8') as bed_input:
                line = bed_input.readline()
                while line:
                    input_gene = line.split('\t')[6].strip()
                    if input_gene == gene:
                        bed_d[gene] += [line.strip()]
                    line = bed_input.readline()

        print('Mitochondrial bed file dictionary created.')
        return bed_d, 'Mitochondrial'

    def output_bed(self, bed_d, MT_bed_d, genelist_fp, bed_abbr):

        dir = os.path.join(script_path,'database','gene_bed',bed_abbr)
        for d in [MT_bed_d,bed_d]:
            for gene in d:
                with open(f'{dir}/{gene}.bed','w',encoding='utf-8') as bed_output:
                    bed_output.write('\n'.join(d[gene]))
        print(f'{bed_abbr} bed files created.')

    



if __name__ == "__main__":
    
    script_path = os.path.dirname(os.path.abspath(__file__))
    
    args = get_parser()

    os.makedirs(os.path.join(script_path,'database','gene_bed','MGIcapture'),exist_ok=True)
    os.makedirs(os.path.join(script_path,'database','gene_bed','cdsplice5'),exist_ok=True)

    genelist = get_genelist()

    create_bed = create_bed()
    update_d, transcript_d, ncbi_d = create_bed.prepare_dict()
    cdsplice5_bed_d = create_bed.cdsplice5(genelist)
    MGIcapture_bed_d = create_bed.MGIcapture(genelist)
    MT_bed_d = create_bed.mitochondria(genelist)
    create_bed.output_bed(cdsplice5_bed_d, MT_bed_d, args.genelist, 'cdsplice5')
    create_bed.output_bed(MGIcapture_bed_d, MT_bed_d, args.genelist, 'MGIcapture')

