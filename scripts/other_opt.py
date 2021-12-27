import numpy as np
import pandas as pd
import numpy.ma as ma
from scipy.spatial import distance_matrix
#from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import umap
REALMIN = np.finfo(float).tiny
def preprocess(datapath, outpath):
    '''
    This function import raw data and fetch the LR pairs appeared in TPM.
    
    '''
    lr_path = datapath + 'LR_pairs.txt'
    lr_df = pd.read_csv(lr_path, sep='\t', delimiter='\t', header=None)
    tpm_path = datapath + 'TPM.txt'
    data = pd.read_csv(tpm_path, sep='\t', delimiter='\t', header=0, index_col=0)
    
    return lr_df, data
    
def denoise(originmat, k):
    '''
    Denoise the affinity matrix, keep only top k cells around each cell.
    '''
    if len(originmat) <= k:
        result = originmat.copy()
    else:
        result = originmat.copy()
        n = len(result)
        for i in range(n):
            row = result[i]
            # descending order
            idx = row.argsort()[::-1]
            # after k, reset to one
            for j in range(k,n):
                result[i,idx[j]] = 0
        result = (result + result.T) / 2
    return result

def csomap_opt(P,num_dims,condition):
    REALMIN = np.finfo(float).tiny
    n = P.shape[0]
    MOMENTUM = 0.5                                 
    FINAL_MOM = 0.8                       
    MOM_SWITCH_ITER = 250
    MAX_ITER = 1001
    Eta = 1000
    MIN_GAIN = 0.01
    
    # Make sure P are set properly 
    if np.any((P < 0)):
        MIN = np.min(P)
        P = P - MIN
        print('P<0')
    else:
        pass
    
    np.fill_diagonal(P, 0)
    P = 0.5 * (P + P.T)
    P = np.maximum(P / np.sum(P), REALMIN)
    const = np.sum(P * np.log(P))
    # If coordinates not given, random initialize. 
    if isinstance(num_dims,int):
        np.random.seed(np.random.randint(9999))
        coord = (np.random.rand(n, num_dims) - 0.5 ) * 50
    else:
        print('Coords are given.')
        coord = num_dims
        R = np.max(abs(coord))
        coord = coord * 50/R  
    
    incs = np.zeros((coord.shape))
    gains = np.ones((coord.shape))
    # run
    for i in range(1,MAX_ITER):
        sum_current = np.sum(coord**2, 1)
        d = sum_current + (sum_current.T + (-2 * np.dot(coord,coord.T))).T
        num = 1 / (1 + d)
        np.fill_diagonal(num,0)
        Q = np.maximum(num / np.sum(num), REALMIN)
        # gradient #
        P_Q = P - Q
        P_Q[(P_Q > 0) & (d <= 0.01)] = -0.01
        L = P_Q * num
        grads = 4 * np.dot((np.diag(np.sum(L, 0)) - L ),coord)

        gains = np.where(np.sign(grads) == np.sign(incs), gains * 0.8, gains)
        gains = np.where(np.sign(grads) != np.sign(incs), gains + 0.2, gains)
        gains[gains < MIN_GAIN] = MIN_GAIN

        incs = MOMENTUM * incs - Eta * (gains * grads)
        coord = coord + incs
        coord = coord - np.mean(coord,0)

        if i == MOM_SWITCH_ITER:
            MOMENTUM = FINAL_MOM
        
        if not (i%100) :
            cost = const - np.sum(P * np.log(Q))
            print('Iteration ',i,': cost is ',cost)
         
        # rescale to 50 #
        R = np.max(abs(coord))
        if condition in 'tight':
            if R > 50 and not (i % 10):
                coord = coord * 50/R
        else:
            if R > 50 and not (i % 1000):
                coord = coord * 50/R  
                
    return cost,coord
    
  
def embedding(aff, out_dir, left_range = 0, right_range = 30, rep=30, dim=3):
    aff = np.array(aff, dtype = 'f')
    mask1 = (aff < 9e-300) & (aff >= 0)
    aff[mask1]=0.1
    np.fill_diagonal(aff,0)
    mask = aff != 0
    aff[mask] = 1 /aff[mask]
    #D = csr_matrix(aff) too less neighbor will occur
    del mask
    for i in range(int(left_range),int(right_range)):
        for j in range(rep):
            pd.DataFrame(umap.UMAP(n_components=dim, metric = "precomputed", n_neighbors=(i+1)*15, random_state = 100*j+3).fit_transform(aff)).to_csv(out_dir + str(i) + '_' + str(j) + '.csv',index = False, header= False, sep = ',')
   
def csomap_sort_cost(A, dim, rep, out_dir):
    for i in range(0,rep):
        cost,coord = csomap_opt(A,dim,'tight')
        pd.DataFrame(coord).to_csv(out_dir +'tsne_dim_'+ str(dim) + '_Rep_' + str(i) + '.csv',index = False, header= False, sep = ',')

def embedding(aff, out_dir, left_range = 0, right_range = 30, rep=30, dim=3):
    aff = np.array(aff, dtype = 'f')
    mask1 = (aff < 9e-300) & (aff >= 0)
    aff[mask1]=0.1
    np.fill_diagonal(aff,0)
    mask = aff != 0
    aff[mask] = 1 /aff[mask]
    #D = csr_matrix(aff) too less neighbor will occur
    del mask
    for i in range(int(left_range),int(right_range)):
        for j in range(rep):
            pd.DataFrame(umap.UMAP(n_components=dim, metric = "precomputed", n_neighbors=(i+1)*15, random_state = 100*j+3).fit_transform(aff)).to_csv(out_dir + str(i) + '_' + str(j) + '.csv',index = False, header= False, sep = ',')
    
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

    
def sparsify_tsne(Q, N, out_dir, rep = 10, percent = 0.4, dim = 3):
    min_cost = 99
    n = len(N)
    method_const = n * percent / np.log(n)
    H, real_per = sparsify(Q, N, percent, method_const)
    rate = real_per / percent
    const = method_const / rate
    print("old const:",method_const, 'new_const:',const)

    for i in range(rep):
        H, real_per = sparsify(Q, N, percent, const)
        for j in range(0,30):
            cost,coord = csomap_opt(H,dim,'tight')
            pd.DataFrame(coord).to_csv(out_dir + 'c_' + str(const) + '_p_' + str(percent) + '_rep_' + str(i) +'_'+ str(j) +'_tsne.csv',index = False, header= False, sep = ',')

    
