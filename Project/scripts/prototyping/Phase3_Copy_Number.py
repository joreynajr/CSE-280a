import os
import pandas as pd 
import math
import numpy as np

# If bed_start and bed_end are -1 that means that the alignment does not map to the interval. 
# Another way you can come to this conclusion is by calculating the fraction of the read that
# overlaps and setting some cutoffs for what you consider VNTR, junction and Non-VNTR. 

# SETTING case specific variables 
case = 1
seq_len = 10000
vntr_len = 65
coverage = 30

# SETTING project specific variables 
proj_dir = '/frazer01/home/joreyna/repos/CSE-280a/Project/'
sample = 'sequence_{}_case_{}'.format(seq_len, case)
coverage_text = 'coverage_{}'.format(coverage)
out_dir = os.path.join(proj_dir, 'output/pipeline/sample/',  sample + '/', coverage_text + '/')











# PARSE the bedtools intersect -wao result 
wao_fn = os.path.join(out_dir, sample + '.wao')   
data = pd.read_table(wao_fn, header=None)
data.columns= ['Map_Template',
'Map_Start',
'Map_End',
'Read_Name',                
'4',
'Read_Orientation',
'6',
'7',
'8',
'9',                
'Map_BPs_Aligned',
'11',
'Bed_Template',
'Bed_Start',
'Bed_End',
'Overlap']
data = data[['Map_Template',
'Map_Start',
'Map_End',
'Read_Name',                
'Read_Orientation',
'Map_BPs_Aligned',
'Bed_Template',                 
'Bed_Start',
'Bed_End',
'Overlap']]
data['Map_Start'] = pd.to_numeric(data['Map_Start'])
data['Map_End'] = pd.to_numeric(data['Map_End'])
data['Map_BPs_Aligned'] = [int(x.replace(',', '')) for x in data['Map_BPs_Aligned']]
data['Bed_Start'] = pd.to_numeric(data['Bed_Start'])
data['Bed_End'] = pd.to_numeric(data['Bed_End'])
data['Overlap'] = pd.to_numeric(data['Overlap'])

# FILTER out reads whose mapping length is not 100 bp's 
data = data[(data['Map_BPs_Aligned'] <= 150) & (data['Map_BPs_Aligned'] >= 50)] 

# ADD a copy number column from the template sequence 
data['Map_Copy_Number'] = data['Map_Template'].str.extract('seq([0-9]*)', expand=True)
data = data[['Map_Copy_Number',
 'Map_Template',
 'Map_Start',
 'Map_End',
 'Read_Name',
 'Read_Orientation',
 'Map_BPs_Aligned',
 'Bed_Template',
 'Bed_Start',
 'Bed_End',
 'Overlap',
 ]]
data = data[data['Overlap'] > 0]

# DETERMINE the copy number 
VNTR_read_count = []
for index, sr in data.iterrows():
    if sr['Overlap'] == sr['Bed_End'] - sr['Bed_Start']:
        VNTR_read_count.append(1)
    else:
        VNTR_read_count.append(float(sr['Overlap']) / (sr['Bed_End'] - sr['Bed_Start']))
data['VNTR_Read_Count'] =  VNTR_read_count
data['Non-VNTR_Read_Count'] = [1 - x for x in VNTR_read_count]
copy_number = data['VNTR_Read_Count'].sum() / data['Non-VNTR_Read_Count'].sum()
print('Copy Number: {}'.format(copy_number))











