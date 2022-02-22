import numpy as np
import pandas as pd
from loess.loess_1d import loess_1d
import lshashpy3 as lshash

REALMIN = np.finfo(float).tiny

def ceil(n):
    return int(-1 * n // 1 * -1)

def floor(n):
    return int(n // 1)

def pear(D,D_re):
    tmp = np.corrcoef(D.flatten(order='C'), D_re.flatten(order='C'))
    return tmp[0,1]  
    
def randomize(mat):
    [m,n] = mat.shape
    # m - spot number
    # n - cell type number
    for i in range(m):
        for j in range(n):
            # loop through each entry
            tmp = mat.iloc[i,j]
            if tmp!=0:
                #print(tmp)
                c = floor(tmp)
                # if entry is integer, pass
                if c == tmp:
                    continue
                else:  
                    d = ceil(tmp)
                    #print(c,d)
                    new = np.random.choice([c,d], p=[d-tmp,tmp-c])
                    mat.iloc[i,j] = new

def data_clean(sc_exp, st_exp):
    # cell x genes
    # 1. remove unexpressed genes
    filtered_sc = sc_exp.loc[:,(sc_exp != 0).any(axis=0)]
    filtered_st = st_exp.loc[:,(st_exp != 0).any(axis=0)]
    st_gene = set(filtered_st.columns)
    sc_gene = set(filtered_sc.columns)
    shared_genes = st_gene.intersection(sc_gene)
    filtered_sc1 = filtered_sc.loc[:,shared_genes]
    filtered_st1 = filtered_st.loc[:,shared_genes]
    return filtered_sc1, filtered_st1 

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
    print('- Given the user-defined parameter keep_lr_per, STORM kept %s, that is, %d of highly variable LR genes.'%(p,take_num))
    a = lr_var.sort_values(by = 0, ascending=False).iloc[0:take_num].index.tolist()
    # combine with top k feature genes
    feature_genes = sort_reg_var.iloc[0:k].index.tolist()
    lr_feature_genes = list(set(feature_genes + a))
    print('- STORM selects %d feature genes.'%(len(lr_feature_genes)))
    return lr_feature_genes

def norm_center(data):
    #first sum to one, then centered
    df = pd.DataFrame(data)
    a = df.apply(lambda x: (x)/np.sum(x) , axis=1)
    return a.apply(lambda x: (x - np.mean(x)) , axis=1)

def init_solution(new_st_exp, new_sc_exp, st_tp_num, meta_df):
# initial solution
    picked_index = {}
    correlations = []
    spot_num = new_st_exp.shape[0]
    gene_num = new_st_exp.shape[1]

    for spot in range(spot_num):
        # maximize correlation
        G = np.array(new_st_exp.iloc[spot]).reshape(1,gene_num)
        cor_st_sc = np.dot(G,new_sc_exp.T)
        
        cor_st_sc_tp = pd.concat([pd.DataFrame(cor_st_sc.T).reset_index(drop=True), meta_df.reset_index(drop=True)], axis=1)
        # redundant change !!
        cor_st_sc_tp.columns = ['cor', 'tp']
        col = list(cor_st_sc_tp.columns)
        B = cor_st_sc_tp.sort_values(by = col[::-1],ascending = False)
        #B = pd.DataFrame(B, index = exp.index )
        a = st_tp_num.iloc[spot]
        b = pd.DataFrame(a[a!=0])

        picked_index[spot] = []
        
        for index, row in b.iterrows():
            #print(row)
            tmp = B[B.iloc[:,1] == index].iloc[0:int(row.iloc[0])].index.tolist()
            picked_index[spot].extend(tmp)
        #break
        correlations.append(pear(G.T,np.sum(new_sc_exp.iloc[picked_index[spot]],axis= 0).values))

    index_cell_dict = {}
    for j in range(new_sc_exp.shape[0]):
        index_cell_dict[j] = new_sc_exp.index[j]
        
    spot_cell_dict = {}
    for group in picked_index:
        spot_cell_dict[group] = []
        for cell in picked_index[group]:
            spot_cell_dict[group].append(index_cell_dict[cell])
            
    return spot_cell_dict, correlations

def replace_cell(sc_exp, st_exp, meta_df, tp_lsh, all_lsh, init_spot_cell_dict, num_nerighbor, mode = 'strict', rand = 3):
# i stand for spot index
    spot_cor = []
    new_picked_index = {}
    for i in range(st_exp.shape[0]):
    #for i in range(30):
        new_picked_index[i] = []
        cell_id = init_spot_cell_dict[i].copy()
        
        G = np.array(st_exp.iloc[i])
        max_cor_rep = 0
        for rep in range(rand):
            cell_id = init_spot_cell_dict[i].copy()
            for j in range(len(cell_id)):
                # replace one cell
                cell = cell_id[j]
                cell_id.remove(cell)
                # sum the other
                C = sc_exp.loc[cell_id]
                sum_C = np.array(np.sum(C,axis = 0))
                # find the difference to be made
                D = G - sum_C

                if mode == 'strict':
                    tp = meta_df.loc[cell][0]
                    nearst = tp_lsh[tp].query(D, num_results = num_nerighbor, distance_func="euclidean")
                    result = len(nearst)
                    candi_group = []
                    if result <= num_nerighbor:
                        # tp number insufficient direct
                        candi_group = meta_df[meta_df['bio_celltype'] == tp].index.tolist()  
                    #print('Its ok to be strict')
                    else:
                        # loop from top neighbors
                        for k in range(result):
                            idx = nearst[k][0][1]
                            candi_group.append(idx)
                else:
                    nearst = all_lsh.query(D, num_results = num_nerighbor, distance_func="euclidean")
                    result = len(nearst)
                    candi_group = []
                    for k in range(result):
                        idx = nearst[k][0][1]
                        candi_group.append(idx)

                candi = sc_exp.loc[candi_group]
                candi_sum = candi + sum_C
                norm_candi_sum = candi_sum.apply(lambda x: (x/np.std(x)) , axis=1)
                candi_cor_list = np.dot(G/np.std(G),norm_candi_sum.T)/len(G)
                candi_max_cor_idx = np.argsort(candi_cor_list)[-1:][0]

                cell_idx = candi_group[candi_max_cor_idx]
                g = sum_C + sc_exp.loc[cell_idx]
                max_cor = np.dot(G/np.std(G),np.array(g/np.std(g)))/len(G)
                #print(i, ":", max_cor)
                tmp_cell_id = cell_id.copy()
                tmp_cell_id.insert(0,cell_idx)     
                cell_id = tmp_cell_id

                if max_cor > max_cor_rep:
                    max_cor_rep = max_cor
                    final_cell_id = cell_id.copy()

        spot_cor.append(max_cor_rep)
        new_picked_index[i].extend(final_cell_id)
    return spot_cor, new_picked_index

def sc_agg(weight_orig, st_exp, meta_df, exp, lr_df, num_per_spot, input_mode, path):
    # copy from original mat
    weight = weight_orig.reset_index()
    del weight[weight.columns[0]]
    # eliminating smaller num
    weight[weight<0.001]=0
    # estimated cell number per spot (can be fractional)
    num = weight * num_per_spot
    # randomize to obtain integer cell-type number per spot
    randomize(num)
    num.to_csv(path + 'cell_type_num_per_spot.csv', index = True, header= True, sep = ',')
    # Filtering
    # 1. keep shared and expressed genes of sc and st data
    print('- %d genes in spatial data, %d genes in single-cell data.'%(st_exp.shape[1],exp.shape[1]))
    filtered_sc, filtered_st = data_clean(exp, st_exp)
    print('- %d shared and expressed genes has been kept.'%(filtered_sc.shape[1]))
    # 2. select feature genes
    sort_genes = feature_sort(filtered_sc, degree = 2, span = 0.3)
    lr_feature_genes = lr_shared_top_k_gene(sort_genes, lr_df, k = 3000, keep_lr_per = 0.8)
    # normalization
    new_sc_exp = norm_center(filtered_sc.loc[:,lr_feature_genes]) 
    #new_st_exp = pd.DataFrame(np.array(norm_center(filtered_st.loc[:,lr_feature_genes])))
    new_st_exp = norm_center(filtered_st.loc[:,lr_feature_genes])
    # initialization
    spot_cell_dict, correlation = init_solution(new_st_exp,new_sc_exp,num,meta_df)
    print('initial solution:',min(correlation), np.mean(correlation), max(correlation))

    num_nerighbor = 200

    all_lsh = lshash.LSHash(1,new_sc_exp.shape[1])
    index = new_sc_exp.index.tolist()
    for j in range(new_sc_exp.shape[0]):
        all_lsh.index(np.array(new_sc_exp.iloc[j]),extra_data=index[j])

    tp_lsh = {}
    for tp in set(meta_df.iloc[:,0]):
        tmp_exp = new_sc_exp[meta_df['bio_celltype'] == tp]
        lsh = lshash.LSHash(1,tmp_exp.shape[1])
        for j in range(tmp_exp.shape[0]):
            lsh.index(np.array(tmp_exp.iloc[j]),extra_data = tmp_exp.index.tolist()[j])
        tp_lsh[tp] = lsh

    spot_cor, new_picked_index = replace_cell(new_sc_exp, new_st_exp, meta_df, tp_lsh, all_lsh, spot_cell_dict, num_nerighbor, mode = input_mode, rand = 2)
    print('swapped solution:',min(spot_cor), np.mean(spot_cor), max(spot_cor))
    # new meta
    picked_index_df = pd.DataFrame()
    for key, value in new_picked_index.items():
        tmp = pd.DataFrame(value)
        tmp[1] = key
        picked_index_df = picked_index_df.append(tmp)
    # single-cell aggregate
    sc_imitator = filtered_sc.loc[picked_index_df[0].tolist()]
    sc_imitator.to_csv(path + 'sc_agg.txt',index = True, header= True, sep = '\t')
    picked_index_df.to_csv(path + 'sc_agg_index.txt',index = False, header= False, sep = '\t')
    pd.DataFrame(correlation).to_csv(path + 'init_cor.txt')
    pd.DataFrame(spot_cor).to_csv(path + 'rep_cor.txt') 
    '''
    for idx,row in picked_index_df.iterrows():
        tmp = meta_df.loc[idx]
        tp = tmp['bio_celltype']
        picked_index_df.loc[idx,'tp'] = tp
    picked_index_df.to_csv(path + 'meta_for_eva.csv')
    '''
    return sc_imitator,picked_index_df
