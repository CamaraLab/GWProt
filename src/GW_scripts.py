
"""
This contains various methods I've made so they're not scattered across different notebooks
copied 4/1/2024
"""

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1' #should be first, before ot

import math
import numpy as np
import random
import csv
import read_pdb
from scipy.spatial.distance import *


from cajal import run_gw, qgw, gw_cython # triangle_ineq_parallel, 

import IdInit




def stride_GW_cell(infile, n = np.inf, scaler = lambda x : x):
    # returns a cython GW_cell from the infile
    coords = read_pdb.get_pdb_coords(infile, n=n)
    stri_distr = np.ones(len(coords))/len(coords)
    ipdm = np.vectorize(scaler)(squareform(pdist(coords)))
    return gw_cython.GW_cell(ipdm,  stri_distr)

def get_nearest_mat(dmat):
    outmat = []
    for row in dmat:
        outrow = np.argsort(row)
        outmat.append(outrow)
    return np.stack(outmat)


# 1 get distance matrix

def dist_csv_to_dist_matrix_no_doubles(
    infile: str, #name of csv file 
    pdb_file_list: list[str], #list of all the pdb files we should expect, matrix is output in that order
    doubled_ids: list[str] = [], # list of ids we should remove and exclude, may or may not be in pdb_file_list
    check: bool = True,
    fun = lambda x: x, #applying a function to the data to get distances
    pdb_files = True #whether it assumes the file names are .pdb and thus drops all other lines
    ): 
    
    #print("test")
    #this function also will deal with cases where the csv has distances between a protein and itself
    
    assert len(doubled_ids) == len(set(doubled_ids))
    for id in doubled_ids:
        if id in pdb_file_list:
            pdb_file_list.remove(id)

    entry_dict = {}
    # protein: index in pdb_file_list
    for i in range(len(pdb_file_list)):
        entry_dict[pdb_file_list[i]]=i

    assert len(pdb_file_list) == len(entry_dict.keys())
    assert set(range(len(pdb_file_list))) == set(entry_dict.values())

    dist_mat = np.zeros((len(pdb_file_list), len(pdb_file_list)))


    
    counter = 1
    doubled_prots = set()#for debugging
    with open(infile, "r", newline="") as infile:
        csv_reader = csv.reader(infile, delimiter=",")

        for line in csv_reader:
            if '.pdb' not in line[0] and pdb_files:
                continue            
            if line[0] not in pdb_file_list or line[1] not in pdb_file_list:
                continue    
            if line[0] in doubled_ids or line[1] in doubled_ids:
                continue
            if line[0] == line[1]:
                doubled_prots.add(line[0])
            counter +=1
            if float(line[2]) <= 0:
                dist_mat[entry_dict[line[0]]][entry_dict[line[1]]] = 0

            else:
                
                dist_mat[entry_dict[line[0]]][entry_dict[line[1]]] = max(fun(float(line[2])),0)
                dist_mat[entry_dict[line[1]]][entry_dict[line[0]]] = max(fun(float(line[2])),0)

    files = set(pdb_file_list)       
    #print(files.difference(doubled_prots))
    #print(doubled_prots.difference(files))
#     print(counter)
#     print(len(pdb_file_list))
#     print(len(list(doubled_prots)))
    if check:
        assert 2* counter == len(pdb_file_list)**2 - len(pdb_file_list)
        assert scipy.spatial.distance.is_valid_dm(dist_mat)
    
    #assert len(dist_mat) ==5127
    
    return dist_mat

def get_nearest_mat(dist_mat):
    #takes in the distance matrix, 
    #outputs matrix whose ith row is the list of the indices in order of closeness to i
    # includes that a point is close to itself
    out_mat = []
    for i in range(dist_mat.shape[0]):
        row = dist_mat[i][:]
        out_row = np.argsort(row)
        out_mat.append(out_row)
    return np.stack(out_mat)



def get_overlap(a,b,c,d):
    # returns the length of the intersection of the intervals [a,b] and [c,d]
    assert a <=b and c <= d
    if b<=c or d <= a:
        return 0
    return min(b,d) - max(a,c)

def id_initial_coupling_unif(m,n):
    # returns a mxn np array that is a coupling of the uniform distributions on n and m, 
    # and is close to the identity permutation array
    if m ==n:
        return np.identity(n)*1/n
        
    P = np.zeros((m,n))
    for i in range(m):
        for j in range(n):
            P[i,j] = get_overlap(i/m, (i+1)/m, j/n, (j+1)/n)
    return P
def random_permutation_initial_coupling_unif(m,n, seed = None):
    if seed:
        np.random.seed(seed)
    P = id_initial_coupling_unif(m,n)
    Q = np.random.permutation(P)
    return(Q)

def id_initial_coupling(a,b):
    # takes in two distributions as 1-dim arrays,
    # outputs a coupling close to the identity matrix
    m = len(a)
    n = len(b)
    a_cdf = [0] #cumulative distribution function a_cdf[i] = sum up through i
    for i in range(m):
        a_cdf.append( a[i] + a_cdf[-1])
    #assert a_cdf[-1] == 1
    
    b_cdf = [0] #cumulative distribution function a_cdf[i] = sum up through i
    for i in range(n):
        b_cdf.append( b[i] + b_cdf[-1])
   # assert b_cdf[-1] == 1

    #print(a_cdf) #testing
    #print(b_cdf)#testing
    
    P = np.zeros((m,n),order = 'C')
    for i in range(m):
        for j in range(n):
            P[i,j] = get_overlap(a_cdf[i],a_cdf[i+1], b_cdf[j],b_cdf[j+1])
    return P

def random_permutation_initial_coupling(a,b, seed = None):
    if seed:
        np.random.seed(seed)
    P = id_initial_coupling(a,b)
    Q = np.random.permutation(P)
    return Q
    
def tensor_coupling(a,b):
    return (a[np.newaxis]).T @ b[np.newaxis]

def unif(n):
    return np.ones((n))*1/n

def GW_identity_init(P1,P2):
    #P1, P2 are GW_cells
    # returns the GW distance using the id_initial_coupling
    a = P1.distribution
    b = P2.distribution
    P = id_initial_coupling(a,b) #to do, check (m,n) vs (n,m)
    C = -2 * P1.dmat @ P @ P2.dmat
    return  gw_cython.gw_cython_init_cost( A = P1.dmat, a = a,  c_A = P1.cell_constant, B = P2.dmat, b = b, c_B = P2.cell_constant, C = C)[1]
    
def GW_tensor_init(P1,P2):
    #P1, P2 are GW_cells
    # returns the GW distance using the id_initial_coupling
    a = P1.distribution
    b = P2.distribution
    #print(a.shape)
    #print(b.shape)
    P =  tensor_coupling(a,b) #possibly wrong direction
    #print(P.shape)
    C = -2 * P1.dmat @ P @ P2.dmat
    return  gw_cython.gw_cython_init_cost( A = P1.dmat, a = a,  c_A = P1.cell_constant, B = P2.dmat, b = b, c_B = P2.cell_constant, C = C)[1]
    

def submat(ar, indices):
    x = np.ix_(indices, indices)
    return ar[x]


def GW_from_coords(prot1_coords, prot2_coords, n = np.inf, scaler = lambda x : x):
    # prot1_coords - (n1,3) array or nested list of the CA coords
    # prot2_coords - (n2,3) array or nested list of the CA coords
    # n - downsampling 

    new_indices1=np.linspace(0, len(prot1_coords), num=min(n, len(prot1_coords)), endpoint=False, dtype=np.int_) # samples n of them, evenly spaced
    new_indices2=np.linspace(0, len(prot2_coords), num=min(n, len(prot2_coords)), endpoint=False, dtype=np.int_) # samples n of them, evenly spaced

    downsampled_p1 = prot1_coords[new_indices1, :]
    downsampled_p2 = prot2_coords[new_indices2, :]

    ipdm1 = squareform(pdist(downsampled_p1))
    ipdm2 = squareform(pdist(downsampled_p2))

    
    GW_cell1 = gw_cython.GW_cell( dmat = np.vectorize(scaler)(ipdm1), distribution = unif(min(n, len(prot1_coords))))
    GW_cell2 = gw_cython.GW_cell( dmat = np.vectorize(scaler)(ipdm2), distribution = unif(min(n, len(prot2_coords))))


    return GW_identity_init(GW_cell1, GW_cell2)



def GW_from_ipdms(ipdm1, ipdm2,n=np.inf, scaler = lambda x : x):

    n1 = len(ipdm1)
    if n1 > n:
        ipdm1 = submat(ipdm1, np.linspace(0, n1, num=n, endpoint=False, dtype=np.int_))

    n2 = len(ipdm2)
    if n2 > n:
        ipdm2 = submat(ipdm2, np.linspace(0, n2, num=n, endpoint=False, dtype=np.int_))
    
    ipdm1 = np.vectorize(scaler)(ipdm1)
    ipdm2 = np.vectorize(scaler)(ipdm2)
    
    GW_cell1 = gw_cython.GW_cell( dmat = (ipdm1), distribution = unif(len(ipdm1)))

    GW_cell2 = gw_cython.GW_cell( dmat = (ipdm2), distribution = unif(len(ipdm2)))


    return GW_identity_init(GW_cell1, GW_cell2)

def run_FGW( D1,pI1, D2,pI2, alpha):
    #ipdm1, pI list 1, ipdm2, pI list 2, alpha
    n1 = len(D1)
    n2 = len(D2)
    
    try:
        assert n1 == len(pI1)
        assert n2 == len(pI2)
    except:
        print(D1.shape, D2.shape, len(pI1), len(pI2))
        assert False
    
    a = np.array([np.array([x]) for x in pI1])
    b = np.array(pI2)
    aa = np.broadcast_to(a,(n1,n2))
    bb = np.broadcast_to(b,(n1,n2))
    M = abs(aa-bb)
    G0 = id_initial_coupling_unif(n1,n2)
    
    d = ot.fused_gromov_wasserstein2(M=M, C1=D1, C2=D2, alpha = alpha, p= unif(n1),q=unif(n2), G0 = G0, loss_fun='square_loss')


    return  0.5 * math.sqrt(d)