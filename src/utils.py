import pandas as pd
import numpy as np
from scipy.spatial import distance_matrix
from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix
from scipy import sparse
import random
########## preprocessing ###################
def pear(D,D_re):
    tmp = np.corrcoef(D.flatten(order='C'), D_re.flatten(order='C'))
    return tmp[0,1] 

def check_weight_sum_to_one(matrix):
    # check if the gene sum is
    check = False
    row = matrix.shape[0]
    if np.sum(np.sum(matrix,axis = 1)) == row:
        check = True
    return check

def check_st_coord(st_coord,st_exp):
    if st_coord.shape[1] >2:
        raise ValueError(
            f'Spatial coordinates expected two dimensional, got {st_coord.shape[1]}.')
    else:
        st_coord.columns = ['x','y']
    if st_coord.shape[0] != st_exp.shape[0]:
        raise ValueError(
            f'Accoding to spatial epxression data, spatial coordinates expected {st_exp.shape[0]} spots, got {st_coord.shape[0]}.')
    return st_coord

def check_st_sc_pair(st_exp, sc_exp):
    if st_exp.columns[0].isupper():
        st_species = 'H'
    else:
        st_species = 'M'
    if sc_exp.columns[0].isupper():
        sc_species = 'H'
    else:
        sc_species = 'M'        
    if not st_species == sc_species:
        st_exp.columns = map(lambda x: str(x).upper(), st_exp.columns)
        sc_exp.columns = map(lambda x: str(x).upper(), sc_exp.columns)
        raise ValueError(
            f'The case of the gene in ST and SC expression data is expected to be the same, got st_species as {st_species}, sc_species as {st_species}.')
    return st_exp, sc_exp

def check_decon_type(weight, meta_df, cell_type_key):
    if len(set(weight.columns).intersection(set(meta_df[cell_type_key]))) != len(set(weight.columns)):
        raise ValueError(
            f'Cell type in weight matrix is different from single-cell meta file.')

def check_spots_match(weight, st_exp, st_coord):
    if len(set(weight.index).difference(set(st_exp.index))) != 0:
        raise ValueError(
            f'Spot index in weight matrix is different from ST expression\'s.')
    if len(set(st_coord.index).difference(set(st_exp.index))) != 0:
        raise ValueError(
            f'Spot index in ST coordinates is different from ST expression\'s.')
