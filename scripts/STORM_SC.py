import argparse
import os
import pandas as pd
import numpy as np
import optimizers
REALMIN = np.finfo(float).tiny

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--sc-file', dest='sc_file', required=True, help='Single-cell candidate library of the corresponding ST tissue')
parser.add_argument('-l', '--lr-file', dest='lr_file', required=True, help='Ligand-receptor pair data')
parser.add_argument('-o', '--out-dir', dest='out_dir', required=False, help='Output file path')

args = parser.parse_args()
if args.out_dir is not None:
    path = args.out_dir + '/'
    if not os.path.exists(path):
        os.makedirs(path)
else:
    path = os.getcwd() + '/'

exp = pd.read_csv(args.sc_file,sep='\t',header=0,index_col = 0, error_bad_lines=False)
lr_df = pd.read_csv(args.lr_file,sep='\t',header=None )
print('file loaded')

# The range for neighbor number
left_range = 0
right_range = 30
# The iteration number for each neighbor
steps = 30
# The number for sparsification repeats
rep = 2
# embedding dimension of the reconstruction
dim = 3
# The percentage of edges to be preserved
percent_list = [0.05,0.2,0.4]

affinitymat = optimizers.calculate_affinity_mat(lr_df, exp)
Q = optimizers.pre_cal(affinitymat)  

for per in percent_list:
    optimizers.sparsify_embedding(Q, affinitymat, path, per, left_range, right_range, steps, rep, dim)


