#! /home/joreyna/anaconda2/envs/hla/bin/python 
import argparse 
import os
import time 
import numpy as np
import copy

project_dir = os.path.join(os.getcwd(), '../')
project_dir = os.path.abspath(project_dir)

# PARSING commandline arguments 
parser = argparse.ArgumentParser(description='Generate a DNA sequence containing a VNTR sequence.')

parser.add_argument('len', metavar='seqLen', type=int, \
		help='The length of the sequences.')

parser.add_argument('vntr', metavar='VNTR', type=str, \
		help='The VNTR that will be introduced.', 
		default='GCACGCTGCTGTGTAGTGGAGAAAGGGCAGGCAGCGAGCAAGCGTGTACAAGGTATATACGTGCC')

parser.add_argument('numVNTR', metavar='numVNTR', type=int, \
		help='The number of VNTR copies that will be introduced.')

parser.add_argument('numMuts', metavar='numMuts', type=int, \
		help='The number of mutations per copy.')

parser.add_argument('--mutation_type', metavar='mutType', type=str, \
		choices=['individual_random_mutations', 'group_random_mutations', 'specific_mutations'], \
		default='individual_random_mutations',
		help='Copies of the VNTR can different mutations. Specify ' + \
			'mutation_type to simulate different mutational ' + \
			'events in the VNTR copies.\n' + \
			'Choices:\n' + \
			'individual_random_mutations,\n' + \
			'group_random_mutations, and\n' + \
			'specific_mutations.')

parser.add_argument('--loc', metavar='locus', type=int, \
		help='The location where the snps are inserted.')

parser.add_argument('--outer_pad', action='store_true', \
		help='Adds a padding around the VNTR for visual aid.', default=False)

parser.add_argument('--inner_pad', action='store_true', \
		help='Adds a padding between copies of the VNTR for visual aid.', default=False)

args = parser.parse_args()


# PRINTING commandline argument values 
print('\n')
print('ArgParse Argument Values')
print('--------------------')
print('len: {}'.format(args.len))
print('VNTR: {}'.format(args.vntr))
print('VNTR copies: {}'.format(args.numVNTR))
print('Mutations per VNTR copy: {}'.format(args.numMuts))
print('Mutation Type: {}'.format(args.mutation_type))
print('location: {}'.format(args.loc))
print('outer pad: {}'.format(args.outer_pad))
print('inner pad: {}'.format(args.inner_pad))
print('\n')


# DEFINING functions for generating random 
# sequences with a VNTR insertion
def generate_mutation(base):
    """
    Taking into account the current base, base, return a mutation.
    
    """
    
    if base == 'A':
        return np.random.choice(['C', 'G', 'T'])
    
    elif base == 'C':
        return np.random.choice(['A', 'G', 'T'])
    
    elif base == 'G':
        return np.random.choice(['A', 'C', 'T'])
    
    elif base == 'T':
        return np.random.choice(['A', 'C', 'G'])
    
    else:
        raise Exception('base is not a proper DNA nucleotide (ACGT).')


def introduce_random_mutations(vntr, m):
    """
    Generate a VNTR sequence with random mutations. The mutations will be the same across different copies. 
    
    Params
    ------
    
    - vntr, the DNA copy sequence which is copied. 
    - m, the number of SNP mutations that will be randomly introduced. 
    
    Returns
    -------
    A single copy of the VNTR sequence with m mutations. \
    """
    
    mutation_sites = np.random.choice(range(len(vntr)), m, replace=False)
    m_vntr = []
    for site, nucleotide in enumerate(vntr):
        if site in mutation_sites:
            m_vntr.append(generate_mutation(nucleotide))
        else:
            m_vntr.append(nucleotide)
    return ''.join(m_vntr)


def introduce_specific_mutations(vntr, sites, mutations):
    """
    Generate a VNTR sequence with the specified mutations at the specified sites. 
    
    Params
    ------
    
    - vntr, the DNA copy sequence which is copied. 
    - sites, locus where the SNP mutation will be introduced. 
    - mutations, a list of mutations.
    
    Returns
    -------
    A single copy of the VNTR sequence with mutations at the specified sites. 
    """
    
    if len(sites) != len(mutations):
        raise Exception('The number of sites and mutations do not correspond.')
        
#     print('VNTR: {}'.format(vntr))
#     print('sites: {}'.format(sites))
#     print('mutations: {}'.format(mutations))
#     print('\n')
    
    m_vntr = list(vntr)
    for site, nucleotide in enumerate(m_vntr):
        
#         print('site: {}'.format(site))
        
        if site in sites:
            mut_idx = sites.index(site)
            
#             print('\tmut_idx: {}, mutant: {}, nucleotide: {}'.\
#                   format(mut_idx, mutations[mut_idx], nucleotide))
            
            if nucleotide == mutations[mut_idx]:
                raise Exception('Not a mutation. The current site is {}. The current '.format(site) +                                 'nucleotide is {}. Please use a different nucleotide '.format(nucleotide) +                                 'for this site.')
            else:
                m_vntr[site] =  mutations[mut_idx]
    return ''.join(m_vntr)


def generate_sequence_with_vntr(sequence, loc, vntr):
    nseq = sequence[0:loc]
    nseq += vntr 
    nseq += sequence[loc:]
    return nseq 

# SETTING a default value for the location 
# of the insert size to the middle of the sequence 
loc = args.loc
if loc == None:
	loc = args.len / 2 

# GENERATE the random sequence 
sequence = ''.join(np.random.choice(['A', 'C', 'G', 'T'], size=args.len))



# MUTATE the vntr copies.
vntr = args.vntr
if args.mutation_type == 'individual_random_mutations':
	# Testing incomplete 
	new_vntr = []
	for i in range(args.numVNTR):
		new_vntr.append(introduce_random_mutations(vntr, args.numMuts))

elif args.mutation_type == 'group_random_mutations':
	# Testing incomplete 
	new_vntr = [introduce_random_mutations(vntr, args.numMuts)] * args.numVNTR 
	
elif args.mutation_type == 'specific_mutations': 
	# Deprecated. Coding incomplete.
	new_vntr = introduce_specific_mutations(vntr, [0], ['C'])


# INSERT inner padding between VNTR copies
if args.inner_pad == True:
	new_vntr = ' '.join(new_vntr)
else:
	new_vntr = ''.join(new_vntr)

# INSERT outer padding around the VNTR
if args.outer_pad == True:
    padding = ' ' * 10
    new_vntr = padding + new_vntr + padding

# INSERT the VNTR into the sequence 
n_sequence = generate_sequence_with_vntr(sequence, loc, new_vntr)


print('Processed Variable Values')
print('--------------------------')
print('sequence: {}'.format(sequence))
print('new_vntr: {}'.format(new_vntr))
print('n_sequence: {}'.format(n_sequence))
print('\n')



# # simulation of 150 bp single ended reads 
# /frazer01/home/joreyna/software/art_src_MountRainier_Linux/examples/run_test_examples_illumina.sh
# To run ART use:
#     /frazer01/home/joreyna/software/art_src_MountRainier_Linux/art_illumina or
#     cd into the directory and say ./art_illumina









