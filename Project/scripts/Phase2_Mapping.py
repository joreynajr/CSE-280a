#! /home/joreyna/anaconda2/envs/hla/bin/python 
import argparse 
import os
import sys
import time 
import copy
import subprocess
import math
import pandas as pd 
import numpy as np

parser = argparse.ArgumentParser(description='Map reads to a reference genome.')
parser.add_argument('--sample', metavar='sampleName', type=str, help='Name of the sample')  
parser.add_argument('--coverage', metavar='coverage', type=int, help='Desired coverage for the sample.')  
parser.add_argument('--copy-len', metavar='vntr-sequence-len', type=int, help='The length of single vntr repeat.')  
parser.add_argument('--non-vntr-len', metavar='vntr-sequence-len', type=int, help='The length of the non-vntr sequence.')  
args = parser.parse_args()

project_dir = os.path.join(sys.argv[0], '../../')
project_dir = os.path.abspath(project_dir)
output_dir = os.path.join(project_dir, 'output/', 'pipeline/', 'sample/', args.sample)
subprocess.call('mkdir -p {}'.format(output_dir), shell=True)

def print_step(step):
	print('*************** {} ***************'.format(step))

#{ RUNNING ART 
#print('\n')
#print_step('Step1: Simulating reads')
art_dir = os.path.join(project_dir, 'output/', 'pipeline/', \
	'sample/', args.sample + '/', 'coverage_{}'.format(args.coverage))
subprocess.call('mkdir -p {}'.format(art_dir), shell=True)
seq_fn=os.path.join(output_dir, '{}.fa'.format(args.sample))
ref_fn=os.path.join(output_dir, '{}_reference.fa'.format(args.sample))
fq_prefix = os.path.join(output_dir, args.sample)
fq_fn = fq_prefix + '.fq' 
cmd='/frazer01/home/joreyna/software/art_src_MountRainier_Linux/art_illumina' + \
	' -ss HS25 -sam -i {} -l 150 -c {} -o {}'.format(seq_fn, args.coverage, fq_prefix)
print('cmd: {}'.format(cmd))
subprocess.call(cmd, shell=True)
#print('\n')
#}

#{ RUNNING alignment using BWA 
#print_step('Step2: Alignment')
bam_fn=os.path.join(art_dir, args.sample + '.bam')
cmd='bwa mem -t 4 {} {} | samtools view -h -b - > {}'.format(ref_fn, fq_fn, bam_fn)
print('cmd: {}'.format(cmd))
subprocess.call(cmd, shell=True)
#print('\n')
#}

#{ INTERSECTING alignments with VNTR's in templates using Bedtools  
#print_step('Step3: Intersecting')
bed_fn=os.path.join(output_dir, args.sample + '_reference.bed') 
wao_fn=os.path.join(art_dir, args.sample + '.wao')
cmd='bedtools intersect -wao -bed -a {} -b {} > {}'.format(bam_fn, bed_fn, wao_fn)
print('cmd: {}'.format(cmd))
subprocess.call(cmd, shell=True)
#print('\n')
#}

#{ SORTING the bam file  
srt_bam_fn = os.path.join(art_dir, args.sample + '.sorted.bam')
cmd='sambamba sort {}'.format(bam_fn)
print('cmd: {}'.format(cmd))
subprocess.call(cmd, shell=True)
#print('\n')
#}


#{ CALCULATE the coverage   
vntr_cov_fn = os.path.join(art_dir, args.sample + '.vntr.cov')
vntr_bed_fn = os.path.join(output_dir, args.sample + '_reference.bed') 
cmd='samtools bedcov {} {} > {}'.format(vntr_bed_fn, srt_bam_fn, vntr_cov_fn)
print('cmd: {}'.format(cmd))
subprocess.call(cmd, shell=True)
#print('\n')

non_vntr_cov_fn = os.path.join(art_dir, args.sample + '.non_vntr.cov')
non_vntr_bed_fn = os.path.join(output_dir, args.sample + '_non_vntr_reference.bed') 
cmd='samtools bedcov {} {} > {}'.format(non_vntr_bed_fn, srt_bam_fn, non_vntr_cov_fn)
print('cmd: {}'.format(cmd))
subprocess.call(cmd, shell=True)
#print('\n')
#}


vntr_depth_fn = os.path.join(art_dir, args.sample + '.vntr.depth')
cmd='samtools depth -b {} {} > {}'.format(vntr_bed_fn, srt_bam_fn, vntr_depth_fn)
print('cmd: {}'.format(cmd))
subprocess.call(cmd, shell=True)

non_vntr_depth_fn = os.path.join(art_dir, args.sample + '.non_vntr.depth')
cmd='samtools depth -b {} {} > {}'.format(non_vntr_bed_fn, srt_bam_fn, non_vntr_depth_fn)
print('cmd: {}'.format(cmd))
subprocess.call(cmd, shell=True)



##{ DETERMINE the copy number  
##print('\n')
##print_step('Step 4: Determining the copy number')
## If bed_start and bed_end are -1 that means that the alignment does not map to the interval. 
## Another way you can come to this conclusion is by calculating the fraction of the read that
## overlaps. 
#
### PARSE the bedtools intersect -wao result 
#data = pd.read_table(wao_fn, header=None)
#data.columns= ['Map_Template',
#'Map_Start',
#'Map_End',
#'Read_Name',                
#'4',
#'Read_Orientation',
#'6',
#'7',
#'8',
#'9',                
#'Map_BPs_Aligned',
#'11',
#'Bed_Template',
#'Bed_Start',
#'Bed_End',
#'Overlap']
#data = data[['Map_Template',
#'Map_Start',
#'Map_End',
#'Read_Name',                
#'Read_Orientation',
#'Map_BPs_Aligned',
#'Bed_Template',                 
#'Bed_Start',
#'Bed_End',
#'Overlap']]
#data['Map_Start'] = pd.to_numeric(data['Map_Start'])
#data['Map_End'] = pd.to_numeric(data['Map_End'])
#data['Map_BPs_Aligned'] = [int(x.replace(',', '')) for x in data['Map_BPs_Aligned']]
#data['Bed_Start'] = pd.to_numeric(data['Bed_Start'])
#data['Bed_End'] = pd.to_numeric(data['Bed_End'])
#data['Overlap'] = pd.to_numeric(data['Overlap'])
#
### FILTER out reads whose mapping length is not 100 bp's 
#data = data[(data['Map_BPs_Aligned'] <= 150) & (data['Map_BPs_Aligned'] >= 50)] 
#
### ADD a copy number column from the template sequence 
#data['Map_Copy_Number'] = data['Map_Template'].str.extract('seq([0-9]*)', expand=True)
#data = data[['Map_Copy_Number',
# 'Map_Template',
# 'Map_Start',
# 'Map_End',
# 'Read_Name',
# 'Read_Orientation',
# 'Map_BPs_Aligned',
# 'Bed_Template',
# 'Bed_Start',
# 'Bed_End',
# 'Overlap',
# ]]
#data = data[data['Overlap'] > 0]
#
### DETERMINE the copy number 
#VNTR_read_count = []
#for index, sr in data.iterrows():
#    if sr['Overlap'] == sr['Bed_End'] - sr['Bed_Start']:
#        VNTR_read_count.append(1)
#    else:
#        VNTR_read_count.append(float(sr['Overlap']) / (sr['Bed_End'] - sr['Bed_Start']))
#data['VNTR_Read_Count'] =  VNTR_read_count
#data['Non-VNTR_Read_Count'] = [1 - x for x in VNTR_read_count]
#
#vntr_read_count_sum = data['VNTR_Read_Count'].sum()
#vntr_read_count_lnorm = vntr_read_count_sum / args.copy_len
#
#non_vntr_read_count_sum = data['Non-VNTR_Read_Count'].sum()
#non_vntr_read_count_lnorm = non_vntr_read_count_sum / args.non_vntr_len
#
#copy_number = str(vntr_read_count_lnorm / non_vntr_read_count_lnorm)
##print('Copy Number: {}'.format(copy_number))
##}
#
#cn_fn = os.path.join(art_dir, args.sample + '.copy_number')
#with open(cn_fn, 'w') as f:
#	f.write('vntr_read_count_sum: {} \n'.format(vntr_read_count_sum))
#	f.write('vntr_read_count_lnorm: {} \n'.format(vntr_read_count_lnorm))
#	f.write('non_vntr_read_count_sum: {} \n'.format(non_vntr_read_count_sum))
#	f.write('non_vntr_read_count_lnorm: {} \n'.format(non_vntr_read_count_lnorm))
#	f.write('copy_number: {} '.format(copy_number))
