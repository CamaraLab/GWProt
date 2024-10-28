
"""
This contains various methods I've made so they're not scattered across different notebooks
copied 4/1/2024
"""

import os
#os.environ['OPENBLAS_NUM_THREADS'] = '1' #should be first, before ot


from threadpoolctl import ThreadpoolController
controller = ThreadpoolController()
import math
import numpy as np
import random
import csv
from scipy.spatial.distance import *


from cajal import run_gw, qgw, gw_cython # triangle_ineq_parallel, 






def get_nearest_mat(dist_mat):
    """
    takes in the distance matrix, 
    outputs matrix whose ith row is the list of the indices in order of closeness to i
    includes that a point is close to itself
    """
    out_mat = []
    for i in range(dist_mat.shape[0]):
        row = dist_mat[i][:]
        out_row = np.argsort(row)
        out_mat.append(out_row)
    return np.stack(out_mat)



def get_overlap(a,b,c,d):
    """returns the length of the intersection of the intervals [a,b] and [c,d]"""
    
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


    
def tensor_coupling(a,b):
    return (a[np.newaxis]).T @ b[np.newaxis]

def unif(n):
    return np.ones((n))*1/n

@controller.wrap(limits=1, user_api='blas')
def GW_identity_init(P1,P2, transport_plan = False):
    #P1, P2 are GW_cells
    # returns the GW distance using the id_initial_coupling
    a = P1.distribution
    b = P2.distribution
    P = id_initial_coupling(a,b) #to do, check (m,n) vs (n,m)
    C = -2 * P1.dmat @ P @ P2.dmat

    res = gw_cython.gw_cython_init_cost( A = P1.dmat, a = a,  c_A = P1.cell_constant, B = P2.dmat, b = b, c_B = P2.cell_constant, C = C,max_iters_ot = 1000000000000)
    if transport_plan:
        return res[1], res[0]
    else:
        return res[1]

@controller.wrap(limits=1, user_api='blas')
def GW_tensor_init(P1,P2, transport_plan = False):
    #P1, P2 are GW_cells
    # returns the GW distance using the id_initial_coupling
    a = P1.distribution
    b = P2.distribution
    #print(a.shape)
    #print(b.shape)
    P =  tensor_coupling(a,b) 
    #print(P.shape)
    C = -2 * P1.dmat @ P @ P2.dmat
    
    res = gw_cython.gw_cython_init_cost( A = P1.dmat, a = a,  c_A = P1.cell_constant, B = P2.dmat, b = b, c_B = P2.cell_constant, C = C,max_iters_ot = 10000000000)
    if transport_plan:
        return res[1], res[0]
    else:
        return res[1]

def submat(ar, indices):
    x = np.ix_(indices, indices)
    return ar[x]





def get_pairing(T, threshold0 = 0.5, threshold1 = 0.5):
    #T is the transport plan matrix
    #returns pairs of indices where T[i,j] > threshold0 * weight0 and T[i,j] > threshold1 * weight1
    
    pairs = []

    weights0 = np.sum(T, axis = 1)
    weights1 = np.sum(T, axis = 0)
    for i in range(T.shape[0]):
        for j in range(T.shape[1]):
            if T[i,j] > threshold0 * weights0[i] and T[i,j] > threshold1 * weights1[j]:
                pairs.append((i,j))
    return pairs



