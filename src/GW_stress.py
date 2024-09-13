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




def GW_stress(prot1, prot2, T):

    n1= len(prot1)
    n2 = len(prot2)
    assert T.shape == (n1,n2)


    A = prot1.ipdm
    a = prot1.distribution
    B = prot2.ipdm
    b = prot2.distribution

    stress1 = np.einsum('ik,il->i',T,(np.einsum('ij,ij->ij',A,A) @T) )   + T @ np.einsum('kl,kl->kl',B,B) @b  -(2 * np.einsum('ab,ab->a', A @T @B, T))
    stress2 = np.einsum('kj,lj->j' ,T @ np.einsum('kl,kl->kl',B,B), T) +a.T @ np.einsum('ij,ij->ij',A,A) @T -(2 * np.einsum('ab,ab->b', A @T @B, T))

    return stress1, stress2
   




def FGW_stress(prot1, prot2, diff_mat, alpha, T):
#now with FGW stress
    n1= len(prot1)
    n2 = len(prot2)
    assert T.shape == (n1,n2)
    assert 0 <= alpha <= 1


    A = prot1.ipdm
    a = prot1.distribution
    B = prot2.ipdm
    b = prot2.distribution



    geo_stress1 = alpha * (np.einsum('ik,il->i',T,(np.einsum('ij,ij->ij',A,A) @T) )   + T @ np.einsum('kl,kl->kl',B,B) @b  -(2 * np.einsum('ab,ab->a', A @T @B, T)))
    geo_stress2 = alpha *(np.einsum('kj,lj->j' ,T @ np.einsum('kl,kl->kl',B,B), T) +a.T @ np.einsum('ij,ij->ij',A,A) @T -(2 * np.einsum('ab,ab->b', A @T @B, T)))

    


    fused_stress1 = (1-alpha) * np.einsum('ik,ik->i', M,T)
    fused_stress2 = (1-alpha) * np.einsum('ik,ik->k', M,T)

    stress1 = geo_stress1 + fused_stress1
    stress2 = geo_stress2 + fused_stress2

    return stress1, stress2




def get_eccentricity(prot, p =2):
    #gets the eccentricity of each point with exponent p
    # https://www.math.ucdavis.edu/~saito/data/acha.read.w12/memoli-gromov-dist.pdf 
    # defn 5.3

    ipdm = prot.ipdm
    distr = prot.distribution

    assert p >0
    n = ipdm.shape[0]    

        
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
            

def earth_mover_distance_stress(X,Y, T = None, distr1 = None, distr2 = None):
    """
    Calculates the earth mover distance stresses between point clouds X and Y. 

    To Do
    """


    # D = ot.dist(coords1, coords2)
    # #print(D.shape)
    # a = np.ones(coords1.shape[0])/coords1.shape[0]
    # b = np.ones(coords2.shape[0])/coords2.shape[0]

    # T = ot.emd(a,b, D)
    # stress = np.einsum('ij,ij ->ij', D,T)
    # cost = np.sum(stress)
    # stress1 = np.sum(stress, axis = 1)
    # stress2 = np.sum(stress, axis = 0)
        




















