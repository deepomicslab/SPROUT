import argparse
import os
import pandas as pd
import numpy as np
import umap
from scipy.spatial import distance_matrix
import optimizers
import preprocessing

REALMIN = np.finfo(float).tiny
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--st-file', dest='st_file', required=True, help='Spatial transcriptomics data')
parser.add_argument('-c', '--st-coordinate', dest='st_coord', required=False, help='Spatial coordinates of the spatial transcriptomics data')
parser.add_argument('-v', '--weight-file', dest='w_file', required=True, help='Deconvoluted ST data')
parser.add_argument('-r', '--sc-file', dest='sc_file', required=True, help='Single-cell candidate library of the corresponding ST tissue')
parser.add_argument('-m', '--meta-file', dest='meta_file', required=True, help='Cell-type annotation of the single-cell candidate library')
parser.add_argument('-l', '--lr-file', dest='lr_file', required=True, help='Ligand-receptor pair data')
parser.add_argument('-o', '--out-dir', dest='out_dir', required=False, help='Output file path')
parser.add_argument('-p','--cell-num-per-spot', dest='num_per_spot', default=10, help='Estimated cell number per spot. Default is 10')
parser.add_argument('-a','--selection_mode', dest='mode', default='strict', help='The choice of either gather cells primarily from the same type (strict) or from all cells (wild) in the candidate library')

args = parser.parse_args()
if args.out_dir is not None:
    path = args.out_dir + '/'
else:
    path = os.getcwd() + '/'

weight = pd.read_csv(args.w_file,sep='\t',header=0, index_col = 0)
st_exp = pd.read_csv(args.st_file,sep='\t',header=0,index_col = 0)
meta_df = pd.read_csv(args.meta_file,sep='\t',header=0,index_col = 0).astype({"bio_celltype": str})
exp = pd.read_csv(args.sc_file,sep='\t',header=0,index_col = 0, error_bad_lines=False)
lr_df = pd.read_csv(args.lr_file,sep='\t',header=None )
print('file loaded')

if args.st_coord is not None:
    st_coord = pd.read_csv(args.st_coord,sep='\t',header=0, index_col = 0)
    st_coord.columns = ['x','y']
else:
    y = pd.DataFrame(st_exp.index)[0]
    st_coord = pd.DataFrame(y.str.split('x',1).tolist(), columns = ['x','y']).astype(float).round()
if st_coord.shape[1] >2:
    print('Spatial coordinates have more than two dimensional, please check the input format.')
if st_coord.shape[0] != st_exp.shape[0]:
    print('Spots number of spatial coordinates did not match with the spatial expression data\'s, please check the input.')

if not st_exp.columns[0].isupper() == exp.columns[0].isupper():
    print('The input of ST and SC files may originate from different species, please check the case of the gene name.')
    st_exp.columns = map(lambda x: str(x).upper(), st_exp.columns)
    exp.columns = map(lambda x: str(x).upper(), exp.columns)
if set(weight.columns).intersection(set(meta_df['bio_celltype'])) != len(set(weight.columns)):
    print('Cell type in weight matrix is unequal to single-cell meta file')

def embedding(sparse_A, ans, path, verbose = True, left_range = 0, right_range = 30, steps = 30, dim = 2):
    aff = np.array(sparse_A, dtype = 'f')
    mask1 = (aff < 9e-300) & (aff >= 0)
    aff[mask1]=0.1
    np.fill_diagonal(aff,0)
    mask = aff != 0
    aff[mask] = 1 /aff[mask]
    #D = csr_matrix(aff) too less neighbor will occur
    del mask
    max_shape = 0
    print('start embedding')
    if verbose == True:
    # save all reconstructed result
        for i in range(int(left_range),int(right_range)):
            for j in range(steps):
                coord = umap.UMAP(n_components=dim, metric = "precomputed", n_neighbors=(i+1)*15, random_state = 100*j+3).fit_transform(aff)
                pd.DataFrame(coord).to_csv(path + str(i) + '_' + str(j) + '.csv',index = False, header= False, sep = ',')
                coord = np.array(coord)
                D_re = distance_matrix(coord, coord)
                cor = optimizers.pear(np.array(ans),D_re)
                if cor > max_shape:
                    max_shape = cor
                    best_in_shape = coord
            print(i,':',max_shape)
    else:
    # only output the best reconstructed result
        for i in range(int(left_range),int(right_range)):
            for j in range(steps):
                coord = umap.UMAP(n_components=dim, metric = "precomputed", n_neighbors=(i+1)*15, random_state = 100*j+3).fit_transform(aff)
                coord = np.array(coord)
                D_re = distance_matrix(coord, coord)
                cor = optimizers.pear(np.array(ans),D_re)
                if cor > max_shape:
                    max_shape = cor
                    best_in_shape = coord
            print(i,':',max_shape)
    pd.DataFrame(best_in_shape).to_csv(path + 'coord_best_in_shape.csv',index = False, header= False, sep = ',')
    print('reached a correlation in shape at:', max_shape)
    return best_in_shape

print('Start to select single-cell aggregates.')
print('Cell num per spot is: %d, mode as %s'%(args.num_per_spot,args.mode))
sc_imitator, picked_index_df = preprocessing.sc_agg(weight, st_exp, meta_df, exp, lr_df, args.num_per_spot, args.mode, path)
print('Single-cell aggregates selection completed.')
exp_T = sc_imitator.T
affinitymat = optimizers.calculate_affinity_mat(lr_df, exp_T)
ans, sc_spot_neighbor = optimizers.dig_hole_spatial_info(st_coord, picked_index_df, r = 2)
#pd.DataFrame(ans).to_csv(path + 'hipp_ans.csv')
#pd.DataFrame(sc_spot_neighbor).to_csv(path + 'hipp_sc_neighbor.csv')
sparse_A = affinitymat * np.array(sc_spot_neighbor)
sparse_A[sparse_A==0] = 0.1
np.fill_diagonal(sparse_A,1)

# The range for neighbor number
left_range = 0
right_range = 30
# The iteration number for each neighbor
steps = 30
# The number for sparsification repeats
dim = 2
# The percentage of edges to be preserved
coord = embedding(sparse_A, ans, path, left_range, right_range, steps, dim)
print('Finished!')
