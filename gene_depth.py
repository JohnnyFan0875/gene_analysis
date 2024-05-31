import os
import sys
import pandas as pd
import gzip
import time
import numpy as np
import argparse
import subprocess

def get_parser():

    parser = argparse.ArgumentParser()
    parser.add_argument('--depth', '-d', type=int, default=20, help='base depth for analysis. default:20')
    parser.add_argument('--output','-o',required=True,type=str,help='output directory')
    parser.add_argument('--genelist', type=str, required=True, help='gene list file. Delimiter "\\n"')
    parser.add_argument('--bamlist', type=str, help='bam list file')
    parser.add_argument('--dir_bed','-b',type=str,required=True,help='directory where bed files locate. bed file name should be <gene>.bed')
    parser.add_argument('--mapping_quality',type=str,required=True,choices=['MQ0','MQ10','MQ20','MQ30'],nargs='+')
    args = parser.parse_args()

    args.dir_bed = os.path.abspath(args.dir_bed)
    args.output = os.path.abspath(args.output)

    return args

def get_bam_file():

    bam_fp = os.path.abspath(args.bamlist)

    if not os.path.exists(bam_fp):
        sys.exit('Bamlist file not found. Program stop.')

    d = dict()
    with open(bam_fp,'r',encoding='utf-8') as bam_input:
        data = bam_input.read().strip()
        for bam in data.split('\n'):
            bam_name, bam_fp = bam.strip().split('\t')
            if not os.path.exists(bam_fp):
                sys.exit(bam_fp + ' not found. Program stop.')
            else:
                d[bam_name] = bam_fp

    return d.values(), d

def get_genelist():

    genelist_fp = os.path.abspath(args.genelist.strip())

    if not os.path.exists(genelist_fp):
        sys.exit('genelist file not found. Program stop.')

    with open(genelist_fp,'r',encoding='utf-8') as genelist_input:
        data = genelist_input.read().strip()
        genelist = [i.strip() for i in data.split('\n')]
                
    return genelist_fp, genelist

def check_bed_file():

    for gene in genelist:
        if not os.path.exists(os.path.join(args.dir_bed,gene+'.bed')):
            sys.exit(gene+'.bed file not found. Program stop.')

class run_program():

    def sambamba(self, bed_abbr):

        for bam_name in d_bam:
            for MQ_element in args.mapping_quality:
                os.makedirs(os.path.join(args.output,'gene_depth'),exist_ok=True)
                os.chdir(args.output)
                cmd1 = f'less {genelist_fp} '+'''|awk '{print "sambamba depth base -t 16 -L ''' + args.dir_bed + '''/"$1".bed -c 0 -q20 -F \\"not duplicate and mapping_quality >='''+MQ_element.replace('MQ','')+'''\\" -o gene_depth/''' + bam_name + '''_"$1"_''' + MQ_element + '''.txt ''' + d_bam[bam_name] + '''"}' |sh '''
                subprocess.call(cmd1, shell=True)
            print(f'{bam_name} sambamba analysis completed.')

    def final_pct(self):

        output_d = dict()
        for bam_name in d_bam:
            for MQ_element in args.mapping_quality:
                for gene in genelist:

                    filepath = os.path.join(args.output, 'gene_depth', f'{bam_name}_{gene}_{MQ_element}.txt')
                    with open(filepath,'r',encoding='utf-8') as txt_input:
                        txt_input_data = txt_input.read().strip()
                        m,n = 0,0
                        for line in txt_input_data.split('\n')[1:]:
                            m+=1
                            if int(line.split('\t')[2])>= args.depth:
                                n+=1
                        pct = str(np.round(n/m,4))
                        output_d['_'.join([bam_name,MQ_element,gene])] = pct

        for MQ in args.mapping_quality:

            with open(os.path.join(args.output, output_dirname+'_'+MQ+'.xls'),'w',encoding='utf-8') as file_output:
                header_li = sorted([f'{i}_{MQ}' for i in d_bam])
                file_output.write('Gene\t'+'\t'.join(header_li)+'\n')

                for gene in genelist:
                    fill = [gene]
                    for header in header_li:
                        bam_name = '_'.join(header.split('_')[:-1])
                        fill.append(output_d['_'.join([bam_name,MQ,gene])])
                    file_output.write('\t'.join(fill) + '\n')
        
        print('Final percentage file created.')



if __name__ == "__main__":

    script_path = os.path.dirname(os.path.abspath(__file__))
    start_time = time.time()

    #define term
    args = get_parser()
    bam_li, d_bam = get_bam_file()
    genelist_fp, genelist = get_genelist()
    output_dirname = os.path.basename(args.output)

    #check gene bed file
    check_bed_file()

    #sambamba run
    run_program = run_program()
    run_program.sambamba('cdsplice5')
    run_program.sambamba('MGIcapture')

    #calculate stat
    run_program.final_pct()

    end_time = time.time()
    print('*'*30+'\nProgram finished. Total spent time %.2fs for this program.'%(end_time-start_time))