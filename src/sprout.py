import tasklogger
import os
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix

from . import utils
from . import preprocessing as prep
from . import optimizers

class SPROUT:
    def __init__(self, st_exp, st_coord, weight, 
                 sc_exp, meta_df, lr_df,
                 cell_type_key = 'celltype',
                 verbose=1, save_path = ''):
        """
        """
        self.st_exp = st_exp
        self.st_coord = st_coord
        self.st_coord.columns = ['x','y']
        self.weight = weight
        self.sc_exp = sc_exp
        self.meta_df = meta_df
        self.lr_df = lr_df
        self.cell_type_key = cell_type_key
        self.save_path = save_path
        self.verbose = verbose
        tasklogger.set_level(verbose)
        if self.save_path != '':
            if not os.path.exists(self.save_path):
                os.makedirs(self.save_path)
        
    def _check_params(self):
        """
        """
        utils.check_st_sc_pair(self.st_exp, self.sc_exp)
        utils.check_decon_type(self.weight, self.meta_df, self.cell_type_key)
        utils.check_spots_match(self.weight, self.st_exp, self.st_coord)
        self.st_coord = utils.check_st_coord(self.st_coord,self.st_exp)
        self.st_exp, self.sc_exp = utils.check_st_sc_pair(self.st_exp, self.sc_exp)

    def select_sc(self, num_per_spot = 10, mode = 'strict', max_rep = 1, repeat_penalty = 10):
        """
        num_per_spot: Estimated cell number per spot. Default is 10
        mode: The choice of either gather cells strictly from the same type (strict) or from all candidate cells (wild) in the candidate library. Default is strict.
        max_rep: Repeat time of the swap procedure.
        repeat_penalty: When one cell has been picked for [THIS] many times, its probability to be picked again decreases by half.
        """
        self.MODE = mode
        self.MAX_REP = max_rep
        self.num_per_spot = num_per_spot
        self.T_HALF = repeat_penalty
        self._check_params()
        self.trans_id_idx = pd.DataFrame(list(range(self.sc_exp.shape[0])), index = self.sc_exp.index)
        # add 22/12/15	
        if self.num_per_spot == 0:	
            self.num = self.weight	
        else:	
            self.num = prep.randomization(self.weight,self.num_per_spot)
        tasklogger.log_start("Filtering")
        self.st_exp, self.sc_exp = prep.data_clean(self.st_exp, self.sc_exp)
        tasklogger.log_complete("Filtering")
        
        tasklogger.log_start("Feature selection")
        self.sort_genes = prep.feature_sort(self.sc_exp, degree = 2, span = 0.3)
        self.lr_hvg_genes = prep.lr_shared_top_k_gene(self.sort_genes, self.lr_df, k = 3000, keep_lr_per = 0.8)
        tasklogger.log_complete("Feature selection")
        
        self.hvg_st_exp = self.st_exp.loc[:,self.lr_hvg_genes]
        self.hvg_sc_exp = self.sc_exp.loc[:,self.lr_hvg_genes]
        norm_hvg_st = prep.norm_center(self.hvg_st_exp)
        norm_hvg_sc = prep.norm_center(self.hvg_sc_exp)
        self.csr_st_exp = csr_matrix(norm_hvg_st)
        self.csr_sc_exp = csr_matrix(norm_hvg_sc)
        
        tasklogger.log_start("Initial cell selection")
        self.spot_cell_dict, self.init_cor, self.picked_time = prep.init_solution(self.num, self.st_exp.index.tolist(),
            self.csr_st_exp,self.csr_sc_exp,self.meta_df[self.cell_type_key],self.trans_id_idx,
            self.save_path, self.T_HALF)
        tasklogger.log_complete("Initial cell selection")
        
        tasklogger.log_start("Swap")
        self.new_picked_index, self.swap_cor, self.after_picked_time = prep.swap_solution(self.csr_st_exp,self.csr_sc_exp, self.spot_cell_dict, 
            self.meta_df,self.cell_type_key,self.picked_time,
            self.save_path,self.lr_hvg_genes,self.trans_id_idx,
            self.MAX_REP, self.MODE, self.T_HALF)
        tasklogger.log_complete("Swap")
        ###### Output files ########
        picked_index_df = pd.DataFrame()
        for key, value in self.new_picked_index.items():
            tmp = pd.DataFrame(value)
            tmp[1] = key
            picked_index_df = picked_index_df.append(tmp)
        picked_index_df = picked_index_df.reset_index()
        del picked_index_df['index']
        picked_index_df[self.cell_type_key] = self.meta_df.loc[picked_index_df[0]][self.cell_type_key].values
        ''''''
        picked_index_df.columns = ['sc_id','spot','celltype']
        picked_index_df['spot'] = picked_index_df['spot'].astype('str')
        ''''''
        self.picked_index_df = picked_index_df
        picked_index_df.to_csv(f'{self.save_path}/sc_agg_meta.tsv',sep = '\t',header=True,index=True)
        cor_res = pd.DataFrame(np.array([self.init_cor, self.swap_cor]).T, columns=['init','swapped'], index = self.st_exp.index)
        cor_res.to_csv(f'{self.save_path}/sc_selection_result.tsv',sep = '\t',header=True,index=True)
        return cor_res,picked_index_df   
    
    def affinity_sparse(self, chunk_size = 12):
        new_st_coord, sc_dis_mat, ans = optimizers.sc_adj_cal(self.st_coord, self.picked_index_df, alpha = 0.01, chunk_size = chunk_size)
        # eliminate neg
        adata = self.sc_exp.sub(self.sc_exp.min(axis = 1), axis=0)
        adata = adata.loc[self.picked_index_df['sc_id']]
        tasklogger.log_start("Affinity")
        sparse_A = optimizers.chunk_cal_aff(adata, sc_dis_mat, self.lr_df)
        tasklogger.log_complete("Affinity")
        sparse_A[sparse_A!=0] = sparse_A[sparse_A!=0] - 0.1
        sparse_A = sparse_A + np.ones(sparse_A.shape) * 0.1
        return sparse_A,ans

    def spatial_recon(self, left_range = 0, right_range = 20, steps = 10, dim = 2,max_dist = 1):
        """
        left_range : int, default: 0. The index range for neighbor number, the actual neighbor number is (i+1)*10
        right_range : int, default: 20. The index range for neighbor number, the actual neighbor number is (i+1)*10
        steps : int, default: 10. The number of repetitions for each neighbor parameter.
        dim : int, default: 2. The embedding dimension of the reconstruction.
        max_dist: float, default:1, recommand in range [0.5,1]. The ratio of the distance between the cells to its spot center to the furthest distance.
        """
        sparse_A,ans = self.affinity_sparse()
        coord = optimizers.embedding(sparse_A, ans, self.save_path, verbose = False, left_range = left_range, right_range = right_range, steps = steps, dim = dim)
        self.picked_index_df = optimizers.center_shift_embedding(coord, self.picked_index_df, max_dist = max_dist)
        new_coord = self.picked_index_df[['adj_UMAP1','adj_UMAP2']]
        self.sparse_A = sparse_A
        self.ans = ans
        self.coord = new_coord
        self.picked_index_df.to_csv(f'{self.save_path}/sc_agg_meta.tsv',sep = '\t',header=True,index=True)
        new_coord.to_csv(f'{self.save_path}/spot_raw_coord_best.tsv',sep = '\t',header=True,index=True)
        return new_coord
