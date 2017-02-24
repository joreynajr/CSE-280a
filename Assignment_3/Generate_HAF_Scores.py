# coding: utf-8


import sys
import time 
import subprocess
import pandas as pd
import string 
from collections import defaultdict



def str_variable(variable, name, tabs=0):
    return  '{}{}: {}'.format('\t' * tabs, name, variable)



def alphabet_generator():
    for letter in string.ascii_lowercase:
        yield letter


# Use msms to generate a populations growing under selection constraint. Generate many
# samples (100) for each of 5 different time points since onset of selection, but plan your
# simulations as it is a forward simulation, and may take some time. Use the following
# parameters: N = 104, s = 0.1, θ =ρ = 250, n = 200. Use the -SI option to sample at different 
# times since onset of selection. Your selection of sampling times should include time before 
# fixation of the favored allele, and times after fixation. 

# Plot the average of HAF scores
# separately for carriers of the favored allele, and the non-carriers as a function of time since
# onset of selection. The x-axis of your plot should be the time since onset of selection, while
# the y-axis is the HAF-score in units of θn. What do you observe?

# ### Run MSMS

def generate_dataset(out_fn='sample.txt', population_size = 10000, s = 0.1, time = 1, favored_site = 0.50000, freq = 0.00001):
    generation = float(time) / (4 * population_size)
    sAa = 2*population_size*s
    sAA = 2 * sAa
    #     freq = float(1)/(2*population_size)
    #     freq = '{:.5f}'.format(freq)
    
    cmd = 'java -jar ' 
    cmd += '/frazer01/home/joreyna/shared_drive/CSE-280a/Assignment_2/msms/lib/msms.jar '
    cmd += '-ms 200 100 -r 250 -t 250 -N {} -SAa {} -SAA {} '.format(population_size, sAa, sAA)
    cmd += '-SF {} {} -Sp 0.50000 -Smark -threads 4 '.format(generation, freq) 
    cmd += '1> /frazer01/home/joreyna/repos/CSE-280a/Assignment_3/{}'.format(out_fn) 
    
    print('cmd: {}'.format(cmd))
    return subprocess.check_call(cmd, shell=True)
    


# ### Load MSMS results 

def load_dataset(fn):
    with open(fn, 'rb') as f:
        cmd = f.readline()
        code = f.readline()
        blank = f.readline()
        drow = blank
        dataset = []
        count = 0
        while drow != '':

            if drow == '\n':
                slashes = f.readline()
                segsites = f.readline()
                positions = f.readline()
                drow = positions 

            else:
                positions = drow.strip().replace('positions: ', '')
                positions = [float(x) for x in positions.split()]
                data = [positions]
                drow = f.readline()
                while drow != '\n':
                    drow = [int(x) for x in drow.strip()]
                    data.append(drow)
                    drow = f.readline()
                dataset.append(data)

    return dataset


# ### Extract the carriers and non-carriers

def rename_and_remove_duplicates_column_names(col_names):
    col_dict = defaultdict(alphabet_generator)
    new_col_names = []
    for col in col_names:
        new_col_names.append('{}_{}'.format(col, col_dict[col].next()))
    return new_col_names


# ### Generate HAF dataframe

def calculate_haf_data(snp_df):
    """
    Generate a list of lists with HAF values.
    """
    global index, col, entry, row_data

    
    frequencies = snp_df.sum()
    haf_data = [] 
    for index in snp_df.index.tolist():
        row_data = []
        for col in snp_df.columns.tolist():
            entry = snp_df.ix[index, col]
            
            try:
                if len(entry) > 1:
                    entry = entry.iloc[0]
            except:
                entry = entry
            
            if entry == 1:
                row_data.append(frequencies[col])
            else:
                row_data.append(0)
        haf_data.append(row_data)
    return haf_data


# Running MSMS
human_time = sys.argv[0]
time = sys.argv[1]
freq = sys.argv[2]
selective_allele = 0.50000
msms_fn = 'results/timepoint_{}.txt'.format(time)
generate_dataset(out_fn = msms_fn, population_size = 10000,
                 s = 0.1, time = time, favored_site = selective_allele, 
                 freq = freq)

# Loading the dataset 
dataset = load_dataset(msms_fn)

# Calculating the HAF matrix  
dataset_df  = []
total_carriers, total_non_carriers = (0,0)
for data in dataset:
    data_df = pd.DataFrame(data[1:], columns=data[0]) 
    data_df.columns = rename_and_remove_duplicates_column_names(data_df.columns.tolist())
    #     print('Are there any duplicate columns?')
    #     print(data_df.columns[data_df.columns.duplicated(keep=False)])
    
    selective_allele_idx = data[0].index(selective_allele)
    #     print('Are we correctly pointing at the selective allele?')
    #     print(data_df.columns.tolist()[selective_allele_idx])
    
    data_df.rename(columns={'0.5_a': 'favored'}, inplace=True)
    #     print('Are we correctly changing the name of the column?')
    #     print(data_df.columns.tolist()[selective_allele_idx])
    dataset_df.append(data_df)
dataset_df = pd.concat(dataset_df)
dataset_df.fillna(value=0, inplace=True)

# Calculating the HAF scores  
carriers_df = dataset_df[dataset_df.ix[:, 'favored'] == 1]
non_carriers_df = dataset_df[dataset_df.ix[:, 'favored'] == 0]

carriers_df = dataset_df[dataset_df.ix[:, 'favored'] == 1]
carrier_haf_df = pd.DataFrame(calculate_haf_data(carriers_df))
carriers_haf_scores = carrier_haf.sum(axis=1)
carriers_haf_scores.to_csv('results/carrier_haf_scores_{}.tsv'.format(human_tp), sep='\t')

non_carriers_df = dataset_df[dataset_df.ix[:, 'favored'] == 0]
non_carrier_haf_df = pd.DataFrame(calculate_haf_data(non_carriers_df))
non_carriers_haf_scores = non_carrier_haf_df.sum(axis=1)
non_carriers_haf_scores.to_csv('results/non_carriers_haf_scores_{}.tsv'.format(human_tp), sep='\t')
























