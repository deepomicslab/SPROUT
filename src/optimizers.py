import pandas as pd
import numpy as np
from scipy.spatial import distance_matrix
from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix
from scipy import sparse
import random
import umap
########## preprocessing ###################
def pear(D,D_re):
    tmp = np.corrcoef(D.flatten(order='C'), D_re.flatten(order='C'))
    return tmp[0,1] 
########## affnity #########
def sc_prep(st_coord, picked_sc_meta):
    # broadcast_st_adj_sc
    st_coord = st_coord.loc[picked_sc_meta['spot'].unique()]
    idx_lst = st_coord.index.tolist()
    idx_dict = {k: v for v, k in enumerate(idx_lst)}
    picked_sc_meta['indice'] = picked_sc_meta['spot'].map(idx_dict)
    coord_dict_x = {v: k for v, k in enumerate(list(st_coord['x']))}
    coord_dict_y = {v: k for v, k in enumerate(list(st_coord['y']))}
    picked_sc_meta['x'] = picked_sc_meta['indice'].map(coord_dict_x)
    picked_sc_meta['y'] = picked_sc_meta['indice'].map(coord_dict_y)
    picked_sc_meta = picked_sc_meta.sort_values(by = 'indice')
    # dist calculation
    sc_coord = picked_sc_meta[['x','y']]
    return st_coord, sc_coord

def sc_adj_cal(st_coord, picked_sc_meta, alpha = 0, chunk_size = 12):
    st_coord, sc_coord = sc_prep(st_coord, picked_sc_meta)
    # alpha = 0 for visum data
    all_x = np.sort(list(set(st_coord.iloc[:,0])))
    unit_len = all_x[1] - all_x[0]
    r = 2 * unit_len + alpha
    
    indicator = lil_matrix((len(sc_coord),len(sc_coord)))
    ans = {}
    n_last_row = 0
    for process_i in range(chunk_size):
        X = np.array_split(np.array(sc_coord), chunk_size)[process_i]
        Y = sc_coord
        chunk = distance_matrix(X,Y)
        ans[process_i] = chunk
        neigh = [np.flatnonzero(d < r) for d in chunk]
        # print(process_i, chunk.shape)
        for i in range(len(neigh)):
            #print(i + n_last_row)
            indicator[i + n_last_row, neigh[i]] = 1
            #break
        n_last_row += chunk.shape[0]
        #break
    return st_coord, indicator, ans

def chunk_cal_aff(adata, sc_dis_mat, lr_df):
    genes = list(adata.columns)
    lr_df = lr_df[lr_df[0].isin(genes) & lr_df[1].isin(genes)]
    gene_index =dict(zip(genes, range(len(genes))))
    index = lr_df.replace({0: gene_index, 1:gene_index}).astype(int)
    ligandindex = index[0].reset_index()[0]
    receptorindex = index[1].reset_index()[1]
    scores = index[2].reset_index()[2]
    Atotake = ligandindex
    Btotake = receptorindex
    allscores = scores
    idx_data = csr_matrix(adata).T
    for i in range(len(ligandindex)):
        if ligandindex[i] != receptorindex[i]:
            Atotake = Atotake.append(pd.Series(receptorindex[i]),ignore_index=True)
            Btotake = Btotake.append(pd.Series(ligandindex[i]),ignore_index=True)
            allscores = allscores.append(pd.Series(scores[i]),ignore_index=True)
    A = idx_data[Atotake.tolist()]
    B = idx_data[Btotake.tolist()]
    full_A = np.dot(csr_matrix(np.diag(allscores)), A).T  
    chunk_size = 20
    cells = list(range(adata.shape[0]))
    affinitymat = np.array([[]]).reshape(0,adata.shape[0])
    affinitymat = csr_matrix(affinitymat)
    #s = time.time()

    for process_i in range(chunk_size):
        #a = time.time() 
        cell_chunk = list(np.array_split(cells, chunk_size)[process_i])
        chunk_A = full_A[cell_chunk]
        chunk_aff = np.dot(chunk_A, B)
        chunk_dis_mat = sc_dis_mat[cell_chunk]
        sparse_A = chunk_dis_mat.multiply(chunk_aff)
        #print(chunk_aff.sum())
        affinitymat = sparse.vstack([affinitymat, sparse_A])
        #b = time.time()
        #print(f'{process_i} done, cost {(b - a):.2f}s.')
    return affinitymat

def prep_aff_umap(aff):
    #aff = np.array(aff, dtype = 'f')
    mask1 = (aff < 9e-300) & (aff >= 0)
    aff[mask1]=0.1
    np.fill_diagonal(aff,0)
    mask = aff != 0
    aff[mask] = 1 /aff[mask]
    #D = csr_matrix(aff) too less neighbor will occur
    del mask
    return aff

def get_center_idx(aff, meta, st_coord):
    # random version
    spot_random_center = []
    for spot in st_coord.index.tolist():
        tmp = meta[meta['spot'] == spot]
        seed = random.sample(tmp.index.tolist(), 1)
        spot_random_center.extend(seed)
    return spot_random_center

########## embedding ##########
def embedding(sparse_A, ans, path, verbose = False, left_range = 0, right_range = 20, steps = 10, dim = 2):
    aff = prep_aff_umap(np.array(sparse_A))
    max_shape = 0
    if verbose == True:
    # save all reconstructed result
        for i in range(int(left_range),int(right_range)):
            for j in range(steps):
                coord = umap.UMAP(n_components=dim, metric = "precomputed", n_neighbors=(i+1)*10, random_state = 100*j+3).fit_transform(aff)
                cor = coord_eva(coord, ans,i,chunk_size = 12)
                pd.DataFrame(coord).to_csv(path + str(i) + '_' + str(j) + '.csv',index = False, header= False, sep = ',')
                if cor > max_shape:
                    max_shape = cor
                    best_in_shape = coord
    else:
    # only output the best reconstructed result
        for i in range(int(left_range),int(right_range)):
            for j in range(steps):
                coord = umap.UMAP(n_components=dim, metric = "precomputed", n_neighbors=(i+1)*10, random_state = 100*j+3).fit_transform(aff)
                coord = np.array(coord)
                cor = coord_eva(coord, ans, chunk_size = 12)
                if cor > max_shape:
                    max_shape = cor
                    best_in_shape = coord
            #print(i,':',max_shape)
    pd.DataFrame(best_in_shape,columns = ['x','y']).to_csv(path + 'raw_coord_best.csv',index = False, header= True, sep = ',')
    print('Reached a highest correlation in shape at:', max_shape)
    return best_in_shape


def center_shift_embedding(sc_coord, sc_meta_orig, max_dist):
    ##### tailored version #####
    sc_meta = sc_meta_orig.copy()
    sc_meta[['UMAP1','UMAP2']] = sc_coord
    umap_core = sc_meta.groupby('spot').mean()[['UMAP1','UMAP2']]
    idx_lst = umap_core.index.tolist()
    idx_dict = {k: v for v, k in enumerate(idx_lst)}
    coord_dict_x = {v: k for v, k in enumerate(list(umap_core['UMAP1']))}
    coord_dict_y = {v: k for v, k in enumerate(list(umap_core['UMAP2']))}
    sc_meta['indice'] = sc_meta['spot'].map(idx_dict)
    sc_meta['core1'] = sc_meta['indice'].map(coord_dict_x)
    sc_meta['core2'] = sc_meta['indice'].map(coord_dict_y)
    # calculating the unit length for the gap between two spot
    x_coors = np.sort(list(set(sc_meta['x'])))
    unit_len = x_coors[1] - x_coors[0]
    spot_space = unit_len/2
    # calculating the scale factor
    core_dist = pd.DataFrame(distance_matrix(sc_coord,umap_core))
    core_dist['spot'] = sc_meta.spot
    max_center_dist = pd.DataFrame(np.diag(core_dist.groupby('spot').max()))
    scale_factor = list(spot_space/(max_center_dist[max_center_dist!=0])[0])
    scale_factor_dict = {v: k for v, k in enumerate(scale_factor)}
    sc_meta['scale_factor'] = sc_meta['indice'].map(scale_factor_dict)
    # print(sc_meta['scale_factor'].head(5))
    # center shift and scale
    tmp = sc_meta[['UMAP1','UMAP2']] - sc_meta[['core1','core2']].values
    tmp1 = tmp*(max_dist * sc_meta[['scale_factor','scale_factor']].values)
    sc_meta[['adj_UMAP1','adj_UMAP2']] = tmp1 + sc_meta[['x','y']].values
    return sc_meta

    
def coord_eva(coord, ans, i,chunk_size = 12):
    cor_all = 0
    for process_i in range(chunk_size):
        X = np.array_split(np.array(coord), chunk_size)[process_i]
        Y = coord
        chunk = distance_matrix(X,Y)
        cor = pear(ans[process_i], chunk)
        # print(cor)
        cor_all += cor
    print(f'Shape correlation when neighbor para is {i} reaches: {cor_all/chunk_size}')
    return cor_all/chunk_size