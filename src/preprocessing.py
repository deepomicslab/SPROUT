import numpy as np
import pandas as pd
from loess.loess_1d import loess_1d
import time
import tasklogger
from . import utils
REALMIN = np.finfo(float).tiny
from scipy.sparse import csr_matrix

def id_to_idx(trans_id_idx,cell_id):
    return list(trans_id_idx.loc[cell_id][0])

def randomize(mat_orig):
    [m,n] = mat_orig.shape
    mat = mat_orig.copy()
    # m - spot number
    # n - cell type number
    for i in range(m):
        for j in range(n):
            # loop through each entry
            tmp = mat.iloc[i,j]
            if tmp!=0:
                #print(tmp)
                c = np.floor(tmp)
                # if entry is integer, pass
                if c == tmp:
                    continue
                else:  
                    d = np.ceil(tmp)
                    #print(c,d)
                    new = np.random.choice([c,d], p=[d-tmp,tmp-c])
                    mat.iloc[i,j] = new
        if mat.iloc[i].sum(axis = 0) == 0:
            # at least one cell
            arg_max = mat_orig.iloc[i].argmax()
            print(f'{i} is all zero,arg_max {arg_max} ')
            mat.iloc[i,arg_max] = 1
    return mat

def randomization(weight,num_per_spot):
    weight_threshold = 0.001
    if not utils.check_weight_sum_to_one(weight):
        # not sum as one
        weight = pd.DataFrame(weight).div(np.sum(weight, axis = 1), axis = 0)
# eliminating small num
    weight[weight < weight_threshold] = 0
    # estimated cell number per spot (can be fractional)
    num = weight * num_per_spot
    # randomize to obtain integer cell-type number per spot
    num = randomize(num)
    # num.to_csv(path + 'cell_type_num_per_spot.csv', index = True, header= True, sep = ',')
    return num

def data_clean(st_exp, sc_exp):
    # cell x genes
    # 1. remove unexpressed genes
    filtered_st = st_exp.loc[:,(st_exp != 0).any(axis=0)]
    filtered_sc = sc_exp.loc[:,(sc_exp != 0).any(axis=0)]
    st_gene = set(filtered_st.columns)
    sc_gene = set(filtered_sc.columns)
    shared_genes = list(st_gene.intersection(sc_gene))
    filtered_st1 = filtered_st.loc[:,shared_genes]
    filtered_sc1 = filtered_sc.loc[:,shared_genes]
    return filtered_st1, filtered_sc1, 

def feature_sort(exp, degree = 2, span = 0.3):
    # 1. input cell x gene
    # exp: gene - row, cell - column
    exp = exp.T
    # 2. calculate mean and var for each gene
    var = np.array(np.log10(exp.var(axis=1) + REALMIN))
    mean = np.array(np.log10(exp.mean(axis=1) + REALMIN))
    # 3. fit model 
    xout, yout, wout = loess_1d(mean, var, frac = span, degree = degree, rotate=False)
    # 4. calculate standaridized value
    exp_center = exp.apply(lambda x: x - np.mean(x), axis=1)
    Z = exp_center.div(yout, axis=0)
    # 5. clipp value by sqrt(N)
    upper_bound = np.sqrt(exp.shape[1])
    Z[Z>upper_bound] = upper_bound
    # 6. sort
    reg_var = pd.DataFrame(Z.var(axis=1))
    sort_reg_var = reg_var.sort_values(by = 0, ascending=False)
    return sort_reg_var

def lr_shared_top_k_gene(sort_reg_var, lr_df, k = 3000, keep_lr_per = 0.8):
    # shared lr genes
    genes = sort_reg_var.index.tolist()
    lr_share_genes = list(set(lr_df[0]).union(set(lr_df[1])).intersection(set(genes)))
    # keep top lr genes
    lr_var = sort_reg_var.loc[lr_share_genes]
    take_num = int(len(lr_var) * keep_lr_per)
    p = "{:.0%}".format(keep_lr_per)
    tasklogger.log_info('- Given the user-defined parameter keep_lr_per, STORM kept %s, that is, %d highly variable LR genes.'%(p,take_num))
    a = lr_var.sort_values(by = 0, ascending=False).iloc[0:take_num].index.tolist()
    # combine with top k feature genes
    feature_genes = sort_reg_var.iloc[0:k].index.tolist()
    lr_feature_genes = list(set(feature_genes + a))
    tasklogger.log_info('- STORM selects %d feature genes.'%(len(lr_feature_genes)))
    return lr_feature_genes

def norm_center(data):
    #first sum to one, then centered
    df = pd.DataFrame(data)
    a = df.apply(lambda x: (x)/np.sum(x) , axis=1)
    return a.apply(lambda x: (x - np.mean(x)) , axis=1)

############ init ##############
def half_life_prob(t,T=10):
    '''
    # When one cell has been picked for T times, 
    # its prob to be picked again decreases by half.
    # T default as 10
    '''
    return (1/2)**(t/T)

def init_solution(cell_type_num, spot_idx, csr_st_exp, csr_sc_exp, meta_df, trans_id_idx, save_path, T_HALF):
    spot_i = -1
    picked_index = {}
    correlations = []
    sc_index = np.array(meta_df.index)
    meta_df = np.array(meta_df)
    picked_time = pd.DataFrame(np.zeros(len(sc_index)), index = sc_index)
    pd.DataFrame(columns = ['spot','time_cost','cor']).to_csv(f'{save_path}/prep_init.log',sep = ',',header=False,mode='a',index=True)
    for spot_name in spot_idx:
        last = time.time()
        spot_i += 1
        prob = half_life_prob(t = picked_time[0].values,T = T_HALF)
        Es = csr_st_exp[spot_i]
        cor_st_sc = Es.dot(csr_sc_exp.T).toarray()[0]
        adj_cor = cor_st_sc * prob
        cor_tp = np.vstack((adj_cor, meta_df, sc_index, cor_st_sc)).T
        sort_cor_tp = cor_tp[cor_tp[:, 0].argsort()][::-1]
        w_i = cell_type_num.iloc[spot_i]
        est_type = pd.DataFrame(w_i[w_i != 0])
        picked_index[spot_name] = []
        for index, row in est_type.iterrows():
            selected_idx = np.where(sort_cor_tp[:,1] == index)[0][0:int(row[spot_name])]
            # modify picked time
            selected_cell_id = list(sort_cor_tp[selected_idx][:,2])
            picked_time.loc[selected_cell_id] += 1
            picked_index[spot_name].extend(selected_cell_id)
        candi_idx = id_to_idx(trans_id_idx, picked_index[spot_name]) 
        agg_exp = csr_sc_exp[candi_idx].sum(axis = 0)
        cor = np.corrcoef(Es.toarray(),np.array(agg_exp))[0,1]
        correlations.append(cor)
        tmp = pd.DataFrame(np.array([[spot_name, time.time()-last, cor]]))  
        tmp.to_csv(f'{save_path}/prep_init.log',sep = ',',header=False,mode='a',index=True)
    f = open(f'{save_path}/prep_init.log','a')
    print(f'init solution: max - {np.max(correlations):.2f}, \
    mean - {np.mean(correlations):.2f}, \
    min - {np.min(correlations):.2f}',file=f)
    f.close()
    print(f'\t Init solution: max - {np.max(correlations):.2f}, \
    mean - {np.mean(correlations):.2f}, \
    min - {np.min(correlations):.2f}')
    return picked_index, correlations, picked_time

#################### # Swap ###########################
def swap_solution(csr_st_exp, csr_sc_exp, picked_index, sc_meta, key, picked_time,save_path, lr_hvg_genes, trans_id_idx, MAX_REP, MODE, T_HALF):
    # TODO Multiprocessing
    all_proces_N = 1
    process_n = 0
    # Initiation
    gene_num = len(lr_hvg_genes)
    s_st_exp = csr_st_exp
    s_sc_exp = csr_sc_exp
    spot_cor = []
    new_picked_index = {}
    sc_meta_dict = {}
    for tp in set(sc_meta[key]):
        idx = sc_meta[sc_meta[key] == tp].index.tolist()
        sc_meta_dict[tp] = dict.fromkeys(idx,0)
    mutual_spots = list(picked_index.keys())
    spot_idx = list(np.array_split(mutual_spots, all_proces_N)[process_n])
    pd.DataFrame(columns = ['spot','time_cost','cor']).to_csv(f'{save_path}/prep_swap.log',sep = ',',header=False,mode='a',index=True)
    # Start swapping
    spot_i = -1
    for spot in spot_idx:
        max_cor_rep = 0
        max_cor = 999
        rep = 0
        spot_i += 1
        last = time.time()
        new_picked_index[spot] = []
        Es = s_st_exp[spot_i]
        norm_Es = csr_matrix(Es/np.std(Es.toarray()))
        spot_cell_id = picked_index[spot].copy()
        while abs(max_cor - max_cor_rep) > 0.001 and rep <= MAX_REP:
            rep +=1
            for i in range(len(spot_cell_id)):
                cell_i = spot_cell_id[i]
                spot_cell_id.remove(cell_i)
                # print(spot_cell_id)
                spot_cell_idx = id_to_idx(trans_id_idx,spot_cell_id)
                spot_remain_mat = s_sc_exp[spot_cell_idx]
                spot_sum = np.array(np.sum(spot_remain_mat,axis = 0))
                if MODE == 'strict':
                # TODO other mode havent finish yet!
                    removed_type = sc_meta.loc[cell_i][key]
                    candi_cell_id = list(sc_meta_dict[removed_type].keys())
                candi_idx = id_to_idx(trans_id_idx, candi_cell_id)    
                candi_mat = s_sc_exp[candi_idx]
                candi_sum = candi_mat + spot_sum
                norm_candi_sum = csr_matrix(candi_sum/np.std(candi_sum,axis = 1))
                candi_cor_list = np.dot(norm_Es,norm_candi_sum.T)/gene_num
                ### 
                prob = half_life_prob(picked_time[0].values,T_HALF)
                adj_cor = candi_cor_list.multiply(prob[candi_idx]).toarray()
                candi_max_cor_idx = np.argsort(adj_cor[0])[-1:][0]
                swaped_idx = candi_idx[candi_max_cor_idx]
                swaped_id = candi_cell_id[candi_max_cor_idx]
                ###        
                new_agg = spot_sum + s_sc_exp[swaped_idx]
                max_cor = np.corrcoef(new_agg,Es.toarray())[0][1]
                #print(i, ":", max_cor)
                tmp_cell_id = spot_cell_id.copy()
                # print(f'max_cor is {max_cor}; max_rep is {max_cor_rep}')
                if max_cor > max_cor_rep:
                    max_cor_rep = max_cor
                    #print(f'insert {swaped_id} to {tmp_cell_id}')
                    tmp_cell_id.insert(0,swaped_id) 
                    picked_time.loc[swaped_id] += 1
                    picked_time.loc[cell_i] -= 1
                else:
                    #print(f'insert {cell_i} back to {tmp_cell_id}')
                    tmp_cell_id.insert(0,cell_i)
                    picked_time.loc[cell_i] += 1
                spot_cell_id = tmp_cell_id
        # print(spot,max_cor_rep)
        tmp = pd.DataFrame(np.array([[spot, time.time()-last, max_cor_rep]]))   
        tmp.to_csv(f'{save_path}/prep_swap.log',sep = ',',header=False,mode='a',index=True)
        spot_cor.append(max_cor_rep)
        new_picked_index[spot] = spot_cell_id
    f = open(f'{save_path}/prep_swap.log','a')
    print(f'swapped solution: max - {np.max(spot_cor):.2f}, \
    mean - {np.mean(spot_cor):.2f}, \
    min - {np.min(spot_cor):.2f}',file=f)
    f.close()
    print(f'\t Swapped solution: max - {np.max(spot_cor):.2f}, \
    mean - {np.mean(spot_cor):.2f}, \
    min - {np.min(spot_cor):.2f}')
    return new_picked_index, spot_cor, picked_time
