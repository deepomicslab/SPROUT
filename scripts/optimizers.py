import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix
#from scipy.sparse import csr_matrix
import umap
REALMIN = np.finfo(float).tiny

def calculate_affinity_mat(lr_df, data):
    '''
    This function calculate the affinity matrix from TPM and LR pairs.
    '''
    # fetch the ligands' and receptors' indexes in the TPM matrix 
    # data.shape = gene * cell
    genes = data.index.tolist()
    lr_df = lr_df[lr_df[0].isin(genes) & lr_df[1].isin(genes)]
    # replace Gene ID to the index of each gene in data matrix #
    gene_index =dict(zip(genes, range(len(genes))))
    index = lr_df.replace({0: gene_index, 1:gene_index}).astype(int)

    ligandindex = index[0].reset_index()[0]
    receptorindex = index[1].reset_index()[1]
    scores = index[2].reset_index()[2]
    
    Atotake = ligandindex
    Btotake = receptorindex
    allscores = scores
    idx_data = data.reset_index()
    del idx_data[idx_data.columns[0]]
    
    for i in range(len(ligandindex)):
        if ligandindex[i] != receptorindex[i]:
            Atotake = Atotake.append(pd.Series(receptorindex[i]),ignore_index=True)
            Btotake = Btotake.append(pd.Series(ligandindex[i]),ignore_index=True)
            allscores = allscores.append(pd.Series(scores[i]),ignore_index=True)

    A = idx_data.loc[Atotake.tolist()]
    B = idx_data.loc[Btotake.tolist()]

    affinitymat = np.dot(np.dot(np.diag(allscores), A).T , B)
    
    return affinitymat
    
def pear(D,D_re):
    tmp = np.corrcoef(D.flatten(order='C'), D_re.flatten(order='C'))
    return tmp[0,1]     

def pre_cal(N):
    np.fill_diagonal(N,0)
    Degree = np.diag(np.sum(N,axis = 1))
    #print("a")
    L = Degree - N
    #print("b")
    L_inv = np.linalg.pinv(L, rcond = REALMIN)
    #print(L_inv)
    #print("c")
    n = len(N)
    Q = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            reff = L_inv[i,i] + L_inv[j,j] - L_inv[i,j] - L_inv[j,i]
            Q[i,j] = N[i,j] * reff
    return Q


def sparsify(Q, N , percent=0.4, const = 1):
    n = len(N)
    # calculate constant C/ep**2 for expected edge k.
    print("const",const)
    P = np.minimum(1, const * np.log(n) * Q)
    #print(P)
    kept = 0
    H = np.ones((n,n)) * REALMIN
    time = 0
    while kept <= percent:
        if time > 20:
            print("Can not reach the requied percentage.")
            break
        for i in range(n):
            for j in range(n):
                rand = np.random.random(1)[0]
                p = P[i,j]
                if rand <= p:
                    H[i,j] = N[i,j] / p
        kept = len(H[H>REALMIN]) / n**2
        time +=1
        #print(H)
        print('kept is:',kept)
    return H,kept

def sparsify_embedding(Q, N, out_dir, percent = 0.2, left_range = 0, right_range = 30, steps = 30, rep = 5 ,dim = 3):
    n = len(N)
    method_const = (n * percent) /np.log(n)
    H, real_per = sparsify(Q, N, percent, method_const)
    rate = real_per / percent
    const = method_const / rate
    print("old const:",method_const, 'new_const:',const)
    for r in range(rep):
        H, real_per = sparsify(Q, N, percent, const)
        H_sym = 0.5 * (H + H.T)
        #print('b')
        mask1 = (H_sym < 9e-300) & (H_sym>0)
        H_sym[mask1]=0.1
        np.fill_diagonal(H_sym,0)
        mask = H_sym != 0
        H_sym[mask] = 1 /H_sym[mask]
        #D = csr_matrix(H_sym)
        del mask
        for i in range(int(left_range),int(right_range)):
            for j in range(steps):
                print("n_neighbors:", (i+1)*15)
                coord = umap.UMAP(n_components=dim, metric = "precomputed", n_neighbors=(i+1)*15, random_state = 100*j+3).fit_transform(H_sym)
                pd.DataFrame(coord).to_csv(out_dir +'Rep_' + str(r) + '_P_'+ str(percent) + '_'+ str(i) + '_' + str(j) + '.csv', index = False, header= False, sep = ',')

    
# ST   
def pass_adj_st_to_sc(spatial_dis,picked_sc_meta):
    spot_sc_num = picked_sc_meta.groupby(1).count()
    sc_spot_neighbor = np.empty((0,int(np.sum(spot_sc_num))))
    for i in range(spatial_dis.shape[0]):
        cell_num_i = spot_sc_num.iloc[i,0]
        row_m = np.empty((cell_num_i,0))
        for j in range(spatial_dis.shape[1]):
            cell_num_j = spot_sc_num.iloc[j,0]
            row_m = np.hstack([row_m, np.full((cell_num_i,cell_num_j), spatial_dis[i,j])])
        sc_spot_neighbor = np.vstack([sc_spot_neighbor, row_m])
    return sc_spot_neighbor

def dig_hole_spatial_info(st_coord, picked_sc_meta, r = 2):
    # r stands for within how many unit lengths to consider as the neighbor.
    all_x = np.sort(list(set(st_coord['x'])))
    unit_len = all_x[1] - all_x[0]
    # broadcast nonadjacency to sc from spatial
    spatial_dis_mat = distance_matrix(st_coord,st_coord)
    spatial_dis_mat[spatial_dis_mat > r * unit_len ] = REALMIN
    # indicator
    spatial_dis_mat[spatial_dis_mat != REALMIN ] = 1
    np.fill_diagonal(spatial_dis_mat,1)
    sc_spot_neighbor = pass_adj_st_to_sc(spatial_dis_mat,picked_sc_meta)
    sc_spot_neighbor[sc_spot_neighbor==REALMIN] = 0
    # answer 
    spatial_dis_mat = distance_matrix(st_coord,st_coord)
    # dont need to fill 1 in diag, cause same spot should be very close hence dis = 0
    ans = pass_adj_st_to_sc(spatial_dis_mat,picked_sc_meta)
    return ans, sc_spot_neighbor
