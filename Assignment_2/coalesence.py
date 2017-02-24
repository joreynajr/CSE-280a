import pandas as pd
import seaborn as sns 
from matplotlib import pyplot as plt
import numpy as np
import math
from scipy import stats
from IPython.display import Image
import re
import sys; sys.path.append('/frazer01/home/joreyna/repos/CSE-280a/Assignment_1/')
import phylo 
from scipy.stats import chisqprob
from scipy.misc import comb
import random

# #### Psudeocode
#-------------Data Structures-------------
#
#Nodes:
#    - identifier (str) <br>
#        ----> leaf nodes will be named idv# 
#        ----> inner nodes will be coalescence nodes named coal# 
#    - height in the tree (int) 
#    - children_ids, list of children identifiers (list) (0 = left child, 1 = right child) 
#    - children_muts, list of mutations each child contains 
#    - children_branches, list of branch lengths for each child
#   
#-------------Algorithm-------------
#(1) Initialize the tree dictionary with the n initial nodes and assign "idv" identifiers.
#(2) Assign the first node as previous coalesence node.
#(3) Create a coalesence node. 
#    (a) Simulate a branch length using the exponential distribution.
#    (b) Assign height as the height of the previous coalesence node + the current branch length.
#    (c) Assign the two children randomly and remove them from the node available for coalesence.
#    (d) For each child calculate the branch length as height of the current coalesence node minus the height of the child.
#    (e) For each child assign mutants by sampling from the poisson distribution (number of mutations = theta * branch length)           and enforce the limited sites assumption to ensure sites are only mutated once.
#(4) Keep track of global variables.
#    (a) Update the previous coalesence node to the newly created noalesence node.
#    (b) Update the sample size to sample size - 1
#    (c) Update the population size to population size - growth function
#(5) Iterate steps (2) - (4) until coalesence is complete 
#
#(Note) Depending on the population size at t[i] we may have a situation where population size will equal or be less than the sample size which is not possible. Instead what will happen is we will have rapid coalesence using the remaing samples. Don't know how often this will happen but no doubt it'll have to be handled. 
# #### Implementation

class node():
    def __init__(self, name, height = 0, children_ids = [], children_muts = [], children_branches = []):
        self.name = name
        self.height = height
        self.children_ids = children_ids
        self.children_muts = children_muts
        self.children_branches = children_branches
    
    def message(self, tabs=0, skip_muts=True):
        tabs_str = '\t' * tabs 
        message = '{}name: {}\n'.format(tabs_str, self.name)
        message += '{}height: {}\n'.format(tabs_str, self.height)
        if len(self.children_ids) > 0:
            for cid, cmuts, cbranch in zip(self.children_ids, self.children_muts, self.children_branches):
                if skip_muts == True:
                    message += '{}child_id: {}, child_branch: {}\n'.format(tabs_str, cid, cbranch)
                else:
                    message += '{}child_id: {}, child_muts: {}, child_branch: {}\n'.format(tabs_str, cid, cmuts, cbranch)
        else: 
            if skip_muts == True:
                message += '{}child_id: [], child_branch: []\n'.format(tabs_str)            
            else:
                message += '{}child_id: [], child_muts: [], child_branch: []\n'.format(tabs_str)

        return message


def print_tree(tree, tabs=0):
    tab_str = '\t' * tabs
    for k, v in tree.items():
        print '{}id: {}, {}'.format(tab_str, k, v.message(tabs = tabs))


def print_dict(tree, tabs=0):
    tab_str = '\t' * tabs
    for k, v in tree.items():
        print '{}id: {}, {}'.format(tab_str, k, v)


def print_avail_nodes(tree, tabs=0):
    tab_str = '\t' * tabs
    for k, v in tree.items():
        print '{}id: {}'.format(tab_str, k)


def get_mutations(num_muts, region, prev_mut_sites):
    muts = []
    while len(muts) < num_muts:
        mut_site = np.random.randint(1, region)
        if mut_site not in prev_mut_sites:
            muts.append(mut_site)
    return muts


def coin_flip(p):
    """
    Flip a coin using a random number generator.
    """
    return True if random.random() < p else False



def simulate_snp_tree_w_exp(sample_size, theta, population_size, alpha, region):
    
    skip_muts = True
    # Initializations of the tree and initial values 
    individuals = ['idv{}'.format(x + 1) for x in range(sample_size)]
    tree = {idv: node(idv) for idv in individuals} 
    avail_nodes = {k:k for k in tree.keys()}
    prev_coal_node = tree[individuals[0]]
    prev_mut_sites = []
    generation_sizes = [population_size]
    coal_lengths = [0]
    branch_len = 1 # Rounding

    # Construction of the coalescence tree
    coal_count = 1
    while sample_size > 1:
                
        coal_prob = comb(sample_size, 2) / (2 * population_size)
        
        # Adding a new epoch/branch length
        if coin_flip(coal_prob):
            
            # Creating a new coalescence node
            coal_node = node('coal{}'.format(coal_count))
            # Initializing children information for printing purposes 
            coal_node.children_ids = [0,0]
            coal_node.children_branches = [0,0]
            coal_node.children_muts = [[0],[0]]        

            coal_node.height = prev_coal_node.height + branch_len
            tree.update({'coal{}'.format(coal_count): coal_node})

            # Choosing the left child
            left_child = random.choice(avail_nodes.keys())
            left_child_node = tree[avail_nodes.pop(left_child)]

            # Choosing the right child
            right_child = random.choice(avail_nodes.keys())
            right_child_node = tree[avail_nodes.pop(right_child)]

            # Choosing the number of mutations for each child
            left_mutuations = []

            # Update available nodes to include the newly created coalescent node 
            avail_nodes.update({'coal{}'.format(coal_count): 'coal{}'.format(coal_count) })

            # Add children ids, branch lengths and mutations 
            coal_node.children_ids = [left_child, right_child]
            coal_node.children_branches = [left_child_node.height + branch_len, right_child_node.height + branch_len]

            # Dropping mutations on the left child branch
            num_muts = coal_node.children_branches[0] * theta
            muts_left = get_mutations(num_muts, region, prev_mut_sites)
            prev_mut_sites.append(muts_left)

            # Dropping mutations on the left child branch
            num_muts = coal_node.children_branches[1] * theta
            muts_right = get_mutations(num_muts, region, prev_mut_sites)
            prev_mut_sites.append(muts_right)    
            coal_node.children_muts = [muts_left, muts_right]

            # Bookkeeping global variables needed for exp dist, new coalescent counts, etc 
            prev_coal_node = coal_node
            sample_size -= 1
            coal_count += 1
            coal_lengths.append(branch_len)
            
        branch_len += 1
        population_size = math.exp(-alpha) * population_size 
        generation_sizes.append(population_size)
        
    return tree, generation_sizes, coal_lengths


def traverse_and_get_mutations(tnode, mutations, data, tree):
    """
    Traverse the coalescent tree and fill data with SNP information.
    """
    
    #     print('node: {}'.format(tnode.name))
    
    # Base case when we have hit a leaf
    if len(tnode.children_ids) == 0:
        data[tnode.name] = mutations 
        return 
    
    for c in range(len(tnode.children_ids)):
        child_muts = tnode.children_muts[c]
        child_node = tree[tnode.children_ids[c]]
        mutations.extend(child_muts)
        traverse_and_get_mutations(child_node, mutations, data, tree)
        mutations = mutations[0: -len(child_muts)]
        
    return 


def simulate_snp_array_w_exp(sample_size, tree): #theta, population_size, alpha, region):
    
    individuals = ['idv{}'.format(x + 1) for x in range(sample_size)]

    # Locating the last node denoting coalescence 
    node_ids = tree.keys()
    node_ids = [x for x in node_ids if re.search('coal', x)]
    node_ids.sort(reverse=True)
    tnode = tree[node_ids[0]]

    # Traversing the coalesence tree to extract mutations 
    mutations = []
    data = {}
    traverse_and_get_mutations(tnode, mutations, data, tree)

    # Constructing a sorted list of all the mutation sites 
    all_mutations = set()
    for node_name, muts in data.items():
        for mut in muts:
            all_mutations.add(mut)
    all_mutations = list(all_mutations)
    all_mutations.sort()

    # Making an indexer between column positions and mutations names 
    mut_indexes = {site: i for i, site in enumerate(all_mutations)}

    # Generating the SNP matrix in list form 
    snp_array = [0] * len(data)
    snp_array = [[0] * len(all_mutations) for x in range(len(data))]

    for i, indv in enumerate(individuals):
        for mut in data[indv]:
            snp_array[i][mut_indexes[mut]] = 1 
            
            
    return snp_array, all_mutations


debug = 'hard'

if debug == 'wf':
    sample_size = 3
    population_size = 1000
    mew = 40/(4 * population_size)
    theta = 100
    alpha = 0
    region = 100000
    
elif debug == 'medium':
    sample_size = 8
    theta = 40
    population_size = 10**6
    theta = 40
    alpha = 2
    region = 100000
    
elif debug == 'hard':
    sample_size = 100
    theta = 40
    population_size = 10**6
    mew = 40/(4 * population_size)
    theta = 40
    alpha = 7/(4* population_size)
    region = 100000
    


coal_tree, generation_sizes, coal_lengths = simulate_snp_tree_w_exp(sample_size, theta, population_size, alpha, region)
snp_array, snp_header = simulate_snp_array_w_exp(sample_size, coal_tree)
snp_df = pd.DataFrame(snp_array, columns=snp_header)#, index=individuals); df

fig = plt.figure(figsize=(8,8))
plt.plot(range(len(generation_sizes)), generation_sizes);
fig.savefig('./generation_size.png')

fig = plt.figure(figsize=(8,8))
plt.plot(range(len(coal_lengths)), coal_lengths);
fig.savefig('./coalesence_lengths_size.png')
