import time
import re
import os
import math
import numpy as np
import random
import ot
import statistics
import numpy.typing as npt
import csv
import itertools as it
from scipy.spatial.distance import *
import multiprocessing
import multiprocess
import pandas as pd
import cajal

import Bio.PDB
from typing import Iterator, Iterable, Optional, TypeVar, Generic

# imports used in compute_gw_copy.ipynb
from Bio.SVDSuperimposer import SVDSuperimposer
#import tmtools
#from tmtools.io import get_structure, get_residue_data
#from tmtools import tm_align
import warnings
warnings.filterwarnings('ignore')


import IdInit
import GW_scripts
import read_pdb
import run_fasta36
# copied 5/30/2024

def GW_stress(ipdm1, ipdm2, distr1 = None, distr2 = None, transport_plan = False):
    n1 = ipdm1.shape[0]
    n2 = ipdm2.shape[0]
    
    if distr1 is None:
        distr1 = np.ones((n1))* (1/n1)
    if distr2 is None:
        distr2 = np.ones(n2)* (1/n2) 
        
    T0 = ot.gromov.gromov_wasserstein(C1 = ipdm1, C2 = ipdm2, p=distr1, q=distr2, log = True, G0 = IdInit.id_initial_coupling(distr1, distr2))
    cost0 = T0[1]['gw_dist']
    T = T0[0]
    #print('gw computed ', cost0)

    A = ipdm1
    a = distr1
    B = ipdm2
    b = distr2

    stress1 = np.einsum('ik,il->i',T,(np.einsum('ij,ij->ij',A,A) @T) )   + T @ np.einsum('kl,kl->kl',B,B) @b  -(2 * np.einsum('ab,ab->a', A @T @B, T))
    stress2 = np.einsum('kj,lj->j' ,T @ np.einsum('kl,kl->kl',B,B), T) +a.T @ np.einsum('ij,ij->ij',A,A) @T -(2 * np.einsum('ab,ab->b', A @T @B, T))


    if transport_plan:
        if cost0 <= 0:
            return 0 , stress1, stress2, T
        else:
            return 0.5* math.sqrt(cost0), stress1, stress2, T
    else:  
        if cost0 <= 0:
            return 0 , stress1, stress2
        else:
            return 0.5* math.sqrt(cost0), stress1, stress2

def GW_stress_randomized(ipdm1, ipdm2, N = 10, distr1 = None, distr2 = None, transport_plan = False):
    current_least_cost = np.inf
    current_best_plan = None

    n1 = ipdm1.shape[0]
    n2 = ipdm2.shape[0]
    
    if distr1 is None:
        distr1 = np.ones((n1))* (1/n1)
    if distr2 is None:
        distr2 = np.ones(n2)* (1/n2) 

    T0 = ot.gromov.gromov_wasserstein(C1 = ipdm1, C2 = ipdm2, p=distr1, q=distr2, log = True, G0 = IdInit.id_initial_coupling(distr1, distr2))
    cost0 = T0[1]['gw_dist']
    T = T0[0]
    current_least_cost = cost0
    current_best_plan = T
    
    
    for i in range(N):
        G0 = IdInit.id_initial_coupling(distr1, distr2)
        for j in range(50):
            G0 += 0.5 * GW_scripts.random_permutation_initial_coupling(distr1, distr2)
            G0 -= 0.5* GW_scripts.random_permutation_initial_coupling(distr1, distr2)
            
        T0 = ot.gromov.gromov_wasserstein(C1 = ipdm1, C2 = ipdm2, p=distr1, q=distr2, log = True, G0 = G0)
        cost1 = T0[1]['gw_dist']
        T1 = T0[0]

        if 0<=cost1 < current_least_cost :
            print('updated, iteration',i)
            current_least_cost = cost1
            current_best_plan = T1

    
    
    A = ipdm1
    a = distr1
    B = ipdm2
    b = distr2
    T = current_best_plan

    stress1 = np.einsum('ik,il->i',T,(np.einsum('ij,ij->ij',A,A) @T) )   + T @ np.einsum('kl,kl->kl',B,B) @b  -(2 * np.einsum('ab,ab->a', A @T @B, T))
    stress2 = np.einsum('kj,lj->j' ,T @ np.einsum('kl,kl->kl',B,B), T) +a.T @ np.einsum('ij,ij->ij',A,A) @T -(2 * np.einsum('ab,ab->b', A @T @B, T))

    if transport_plan:
        return 0.5* math.sqrt(current_least_cost), stress1, stress2, T
    else:
        
        return 0.5* math.sqrt(current_least_cost), stress1, stress2


def FGW_stress(ipdm1, ipdm2, diff_mat, alpha, distr1 = None, distr2 = None, transport_plan = False):
#now with FGW stress
    n1 = ipdm1.shape[0]
    n2 = ipdm2.shape[0]
    
    if distr1 is None:
        distr1 = np.ones((n1))* (1/n1)
    if distr2 is None:
        distr2 = np.ones(n2)* (1/n2) 
        
    T0 = ot.gromov.fused_gromov_wasserstein(C1 = ipdm1, C2 = ipdm2, M = diff_mat, alpha = alpha, p=distr1, q=distr2, log = True, G0 = IdInit.id_initial_coupling(distr1, distr2))
    cost0 = T0[1]['fgw_dist']
    T = T0[0]
    #print('fgw computed ', cost0)

    A = ipdm1
    a = distr1
    B = ipdm2
    b = distr2

    geo_stress1 = alpha * (np.einsum('ik,il->i',T,(np.einsum('ij,ij->ij',A,A) @T) )   + T @ np.einsum('kl,kl->kl',B,B) @b  -(2 * np.einsum('ab,ab->a', A @T @B, T)))
    geo_stress2 = alpha *(np.einsum('kj,lj->j' ,T @ np.einsum('kl,kl->kl',B,B), T) +a.T @ np.einsum('ij,ij->ij',A,A) @T -(2 * np.einsum('ab,ab->b', A @T @B, T)))

    


    fused_stress1 = (1-alpha) * np.einsum('ik,ik->i', M,T)
    fused_stress2 = (1-alpha) * np.einsum('ik,ik->k', M,T)

    stress1 = geo_stress1 + fused_stress1
    stress2 = geo_stress2 + fused_stress2



    if transport_plan:
        if cost0 <= 0:
            return 0 , stress1, stress2, T
        else:
            return 0.5* math.sqrt(cost0), stress1, stress2, T
    else:  
        if cost0 <= 0:
            return 0 , stress1, stress2
        else:
            return 0.5* math.sqrt(cost0), stress1, stress2

def FGW_stress_seq_from_prots(p1,p2, alpha, transport_plan = False, n= 200, allow_mismatch = True, redo_calc = False):
    # redo_calc - same settings as original calculation

    if redo_calc:
        n = 200
        alpha = 0.5
        allow_mismatch = True
        p1.scale_ipdm(scaler = math.sqrt, inplace = True)
        p2.scale_ipdm(scaler = math.sqrt, inplace = True)
        p1.convolve_pIs_fasta(kernel_list = [1,2,3,4,3,2,1], origin= 3, inplace=True )
        p2.convolve_pIs_fasta(kernel_list = [1,2,3,4,3,2,1], origin= 3, inplace=True )
        
    inds1, inds2 = FGW_protein.FGW_protein.run_ssearch_indices(p1 = p1, p2 = p2,allow_mismatch = allow_mismatch)

    if n < len(inds1):
        l,s = np.linspace(0, len(inds1), num=n, endpoint=False, dtype=int,retstep = True) 
        subindices = np.array([int(i + s//2) for i in l])
        inds1 = [inds1[i] for i in subindices]
        inds2 = [inds2[i] for i in subindices]

    
    p3 = p1.downsample_by_indices(inds1)
    p4 = p2.downsample_by_indices(inds2)

    D3 = p3.ipdm
    D4 = p4.ipdm
    pI3 = p3.pI_list
    pI4 = p4.pI_list
    n3 = len(D3)
    n4 = len(D4)
    distr1 = GW_scripts.unif(n3)
    distr2 = GW_scripts.unif(n4)
    
    try:
        assert n3 == len(pI3)
        assert n4 == len(pI4)
    except:
        print(D3.shape, D4.shape, len(pI3), len(pI4))
        assert False
    
    a = np.array([np.array([x]) for x in pI3])
    b = np.array(pI4)
    aa = np.broadcast_to(a,(n3,n4))
    bb = np.broadcast_to(b,(n3,n4))
    M = abs(aa-bb)
    G0 = GW_scripts.id_initial_coupling_unif(n3,n4)
    
    #d = ot.fused_gromov_wasserstein2(M=M, C1=D1, C2=D2, alpha = alpha, p= GW_scripts.unif(n1),q=GW_scripts.unif(n2), G0 = G0, loss_fun='square_loss')

    T0 = ot.gromov.fused_gromov_wasserstein(C1 = D3, C2 = D4, M = M, alpha = alpha, p=distr1, q=distr2, log = True, G0 = G0)
    cost0 = T0[1]['fgw_dist']
    T = T0[0]
    #print('fgw computed ', cost0)

    A = D3
    a = distr1
    B = D4
    b = distr2

    geo_stress1 = alpha * (np.einsum('ik,il->i',T,(np.einsum('ij,ij->ij',A,A) @T) )   + T @ np.einsum('kl,kl->kl',B,B) @b  -(2 * np.einsum('ab,ab->a', A @T @B, T)))
    geo_stress2 = alpha *(np.einsum('kj,lj->j' ,T @ np.einsum('kl,kl->kl',B,B), T) +a.T @ np.einsum('ij,ij->ij',A,A) @T -(2 * np.einsum('ab,ab->b', A @T @B, T)))
    


    fused_stress1 = (1-alpha) * np.einsum('ik,ik->i', M,T)
    fused_stress2 = (1-alpha) * np.einsum('ik,ik->k', M,T)

    stress1 = geo_stress1 + fused_stress1
    stress2 = geo_stress2 + fused_stress2

    if transport_plan:
        if cost0 <= 0:
            return 0 , stress1, stress2,inds1, inds2, T
        else:
            return 0.5* math.sqrt(cost0), stress1, stress2, inds1, inds2, T
    else:  
        if cost0 <= 0:
            return 0 , stress1, stress2, inds1, inds2,
        else:
            return 0.5* math.sqrt(cost0), stress1, stress2, inds1, inds2

def FGW_stress_from_prots(p1,p2, alpha, transport_plan = False):
    # redo_calc - same settings as original calculation

    
        
    p3 = p1
    p4 = p2

    D3 = p3.ipdm
    D4 = p4.ipdm
    pI3 = p3.pI_list
    pI4 = p4.pI_list
    n3 = len(D3)
    n4 = len(D4)
    distr1 = GW_scripts.unif(n3)
    distr2 = GW_scripts.unif(n4)
    
    try:
        assert n3 == len(pI3)
        assert n4 == len(pI4)
    except:
        print(D3.shape, D4.shape, len(pI3), len(pI4))
        assert False
    
    a = np.array([np.array([x]) for x in pI3])
    b = np.array(pI4)
    aa = np.broadcast_to(a,(n3,n4))
    bb = np.broadcast_to(b,(n3,n4))
    M = abs(aa-bb)
    G0 = GW_scripts.id_initial_coupling_unif(n3,n4)
    
    #d = ot.fused_gromov_wasserstein2(M=M, C1=D1, C2=D2, alpha = alpha, p= GW_scripts.unif(n1),q=GW_scripts.unif(n2), G0 = G0, loss_fun='square_loss')

    T0 = ot.gromov.fused_gromov_wasserstein(C1 = D3, C2 = D4, M = M, alpha = alpha, p=distr1, q=distr2, log = True, G0 = G0)
    cost0 = T0[1]['fgw_dist']
    T = T0[0]
    #print('fgw computed ', cost0)

    A = D3
    a = distr1
    B = D4
    b = distr2

    geo_stress1 = alpha * (np.einsum('ik,il->i',T,(np.einsum('ij,ij->ij',A,A) @T) )   + T @ np.einsum('kl,kl->kl',B,B) @b  -(2 * np.einsum('ab,ab->a', A @T @B, T)))
    geo_stress2 = alpha *(np.einsum('kj,lj->j' ,T @ np.einsum('kl,kl->kl',B,B), T) +a.T @ np.einsum('ij,ij->ij',A,A) @T -(2 * np.einsum('ab,ab->b', A @T @B, T)))

    


    fused_stress1 = (1-alpha) * np.einsum('ik,ik->i', M,T)
    fused_stress2 = (1-alpha) * np.einsum('ik,ik->k', M,T)

    assert (M >= 0).all()
    assert (T >=0).all()

    assert (fused_stress1 > -1e-1).all() #1e-5 is for margin of floating point error
    assert (fused_stress2 > -1e-1).all()
    
    stress1 = geo_stress1 + fused_stress1
    stress2 = geo_stress2 + fused_stress2

    if transport_plan:
        if cost0 <= 0:
            return 0 , stress1, stress2, T
        else:
            return 0.5* math.sqrt(cost0), stress1, stress2, inds1, inds2, T
    else:  
        if cost0 <= 0:
            return 0 , stress1, stress2, inds1, inds2,
        else:
            return 0.5* math.sqrt(cost0), stress1, stress2, 


def GW_stress_from_prots(p1,p2, transport_plan = False):
    #do not downsample here!
    return GW_stress(p1.ipdm, p2.ipdm, transport_plan = transport_plan)


def GW_stress_from_prots_randomized(p1,p2, N=10, transport_plan = False):
    #do not downsample here!
    return GW_stress_randomized(p1.ipdm, p2.ipdm, transport_plan = transport_plan)

def get_eccentricity(ipdm, p =2, distr = None):
    #gets the eccentricity of each point with exponent p
    # https://www.math.ucdavis.edu/~saito/data/acha.read.w12/memoli-gromov-dist.pdf 
    # defn 5.3
    assert p >0
    n = ipdm.shape[0]    
    if distr is None:
        distr = np.ones((n))* (1/n)
        
    if p == np.inf:
        ipdm_pp =  ipdm* (distr != 0)
        eccentricity = np.max(ipdm_pp, axis = 0)
    else:
        ipdm_pp = ipdm**p
        ipdm_pp_w = ipdm_pp * distr #not sure about shape and broadcasting here
        pre_stress = np.sum(ipdm_pp_w, axis = 0) #unclear which axis this should be
        eccentricity = pre_stress**(1/p)

    
    return eccentricity
    
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
            
        
        
  


















