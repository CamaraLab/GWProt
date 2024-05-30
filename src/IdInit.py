

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1' #should be first, before ot

#copied 4/1/2024

import math
import numpy as np
import random
import read_pdb

from cajal import run_gw, qgw,  gw_cython # triangle_ineq_parallel,



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
    P = id_initial_coupling(m,n)
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
    
    P = np.zeros((m,n))
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
    
def GW_init(P1,P2, init_type, seed = None):
    #P1, P2 are GW_cells
    # returns the GW distance using the id_initial_coupling
    a = P1.distribution
    m = len(a)
    b = P2.distribution
    n = len(b)
    match init_type:
        case 'id':
            P = id_initial_coupling(a,b) #to do, check (m,n) vs (n,m)
        case 'id_rand':
            P = random_permutation_initial_coupling(a,b,seed)
            if (P == id_initial_coupling(a,b)).all():
                print('random permutation is the identity')
        case 'unif':
            P = np.ones((m,n))/(n*m)
        case 'tensor':
            P = tensor_coupling(a,b)
    
    C = -2 * P1.dmat @ P @ P2.dmat
    return  gw_cython.gw_cython_init_cost( A = P1.dmat, a = a,  c_A = P1.cell_constant, B = P2.dmat, b = b, c_B = P2.cell_constant, C = C)[1]
    



#test 
