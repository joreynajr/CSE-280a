from glob import glob
import string 
import copy
from itertools import product

def load_snp_data(fn):
    """Load the snp data into an array (list of lists)"""
    snp_array = []
    with open(fn, 'r') as data:
        for line in data:
            row = []
            line = line.strip()
            for entry in line: 
                row.append(int(entry))
            snp_array.append(row)
    return snp_array

def transpose_array(array):
    array = copy.deepcopy(array)
    t_array = []
    for column in range(len(array[0])):
        t_column = []
        for row in range(len(array)):
            t_column.append(array[row][column])
        t_array.append(t_column)
    return t_array

def get_intersects(column_i, column_j):
    """
    Looking for mutations in column_j that intersect with 
    mutations in column_i.
    """
    intersects = []
    for index in range(len(column_i)):   
        if column_j[index] == 1 and column_i[index] == 1:
            intersects.append(index)
    return intersects


def get_differences(column_i, column_j):
    """
    Looking for mutations in column_j not found in column_i.
    """    
    differences = []
    for index in range(len(column_i)):   
        if column_j[index] == 1 and column_i[index] == 0:
            differences.append(index)
    return differences

def count_mutations(column):
    """
    Counting the number of mutations at a given site (column)."""
    mutations = 0
    for entry in column:
        if entry == 1:
            mutations += 1
    return mutations 


def is_perfect_phylo(t_array):
    
    column = 0 
    perfect_phylo = True
    num_columns = len(t_array)
    num_rows = len(t_array[0])
    while column < num_columns - 1:    


        message = 'column: {}'.format(column)
        column_i = t_array[column]
        column_j= t_array[column + 1]

        mutations = count_mutations(column_j)
        intersects = get_intersects(column_i, column_j)
        differences = get_differences(column_i, column_j)

        if len(intersects) == mutations:
            # column j is a subset. 
            column += 1

        elif len(differences) > 0 and len(intersects) == 0:
            # column j is not a subset and does not share 
            # any mutations with column i.
            column += 1

        elif len(differences) > 0 and len(intersects) > 0:
            # column j is only a partial subset of column i
            # so their intersection is not trivial and the 
            # perfect phylogeny is lost.
            perfect_phylo = False
            break


        message += ', mutations: {}, intersects: {}, differences: {}'.\
            format(mutations, intersects, differences)
        #time.sleep(1)
        #print message
        #print 
        
    return perfect_phylo
