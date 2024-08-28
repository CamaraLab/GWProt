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


import sys
sys.path.insert(0,'../PGC020.a12/src')



import GW_scripts
import read_pdb
import FGW_protein
import IdInit
import pymol_tmalign_wrapper_Copy1
from pymol import cmd



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

def FGW_stress_from_prots(p1,p2, alpha, transport_plan = False, n= np.inf, d = None):
    #d is the aa distance dictionary


    if d is None:
        #use  pI data
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
    else:
        D3 = p1.ipdm
        D4 = p2.ipdm
        n1 = len(p1)
        n2 = len(p2)
        n3 = n1
        n4 = n2
        distr1 = GW_scripts.unif(n3)
        distr2 = GW_scripts.unif(n4)
        # M has shape (n1,n2)
        M = np.zeros((n1,n2))
        for i in range(n1):
            for j in range(n2):
                M[i,j] = d[p1.seq[i]][p2.seq[j]]

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

    # assert (geo_stress1 > -1e-1).all() #this often fails
    # assert (geo_stress2 > -1e-1).all()


    
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
            
        



def get_pymol_transport(file1, file2, aligner = 'cealign', TM_exe = '/home/elijah/pymol/bin/TMalign' ): #filepaths to pbds
    # TM_exe filepath to TM align
    cmd.delete('all')
    cmd.load(file1, 'prot1')
    cmd.load(file2, 'prot2')

    coords01 = []
    coords02 = []
    cmd.iterate_state(1, "prot1 and name CA", "coords01.append((x, y, z))", space={'coords01': coords01})
    cmd.iterate_state(1, "prot2 and name CA", "coords02.append((x, y, z))", space={'coords02': coords02})
    coords01 = np.stack(coords01)
    coords02 = np.stack(coords02)

    match aligner:
        case 'cealign':
            cmd.cealign('prot1', 'prot2')
        case 'align':
            cmd.align('prot1', 'prot2')
        case 'super':
            cmd.super('prot1', 'prot2')     
        case 'tmalign':
            pymol_tmalign_wrapper_Copy1.tmalign('prot1', 'prot2', quiet=1, exe = TM_exe, return_alignment=False) 
            #cmd.tmalign('prot1', 'prot2') 
        case _:
            raise ValueError("valid arguments for aligner are 'align', 'cealign', 'super', tmalign")
    coords1 = []
    coords2 = []
    cmd.iterate_state(1, "prot1 and name CA", "coords1.append((x, y, z))", space={'coords1': coords1})
    cmd.iterate_state(1, "prot2 and name CA", "coords2.append((x, y, z))", space={'coords2': coords2})
    coords1 = np.stack(coords1)
    coords2 = np.stack(coords2)
    
    
    #print(coords1.shape)
    #print(coords2.shape)
    D = ot.dist(coords1, coords2)
    #print(D.shape)
    a = np.ones(coords1.shape[0])/coords1.shape[0]
    b = np.ones(coords2.shape[0])/coords2.shape[0]

    T = ot.emd(a,b, D)
    stress = np.einsum('ij,ij ->ij', D,T)
    cost = np.sum(stress)
    stress1 = np.sum(stress, axis = 1)
    stress2 = np.sum(stress, axis = 0)

    return T, cost, stress1, stress2 


#warning -  the pymol coords could be returned in a different order
# or try cmd.get_coords(str selection)
# not totally sure how ordering works

# pymol_tmalign_wrapper_Copy1.tmalign('prot1', 'prot2', quiet=1, return_alignment=True)
    
    


















