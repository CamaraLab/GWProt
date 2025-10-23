import os
import math
import numpy as np
import ot
import numpy.typing as npt
from scipy.spatial.distance import *
from typing import Iterator, Iterable, Optional, TypeVar, Generic


 


from .GW_scripts import *
from .GW_protein import *




def GW_lgd(prot1:GW_protein, prot2:GW_protein, T: np.array):
    """
    :param prot1:
    
    """

    n1= len(prot1)
    n2 = len(prot2)
    assert T.shape == (n1,n2)


    A = prot1.ipdm
    a = prot1.distribution
    B = prot2.ipdm
    b = prot2.distribution

    lgd1 = np.einsum('ik,il->i',T,(np.einsum('ij,ij->ij',A,A) @T) )   + T @ np.einsum('kl,kl->kl',B,B) @b  -(2 * np.einsum('ab,ab->a', A @T @B, T))
    lgd2 = np.einsum('kj,lj->j' ,T @ np.einsum('kl,kl->kl',B,B), T) +a.T @ np.einsum('ij,ij->ij',A,A) @T -(2 * np.einsum('ab,ab->b', A @T @B, T))

    return lgd1, lgd2
   




def FGW_lgd(prot1, prot2, diff_mat, alpha, T):
#now with FGW lgd
    n1= len(prot1)
    n2 = len(prot2)
    assert T.shape == (n1,n2)
    assert 0 <= alpha <= 1


    A = prot1.ipdm
    a = prot1.distribution
    B = prot2.ipdm
    b = prot2.distribution



    geo_lgd1 = alpha * (np.einsum('ik,il->i',T,(np.einsum('ij,ij->ij',A,A) @T) )   + T @ np.einsum('kl,kl->kl',B,B) @b  -(2 * np.einsum('ab,ab->a', A @T @B, T)))
    geo_lgd2 = alpha *(np.einsum('kj,lj->j' ,T @ np.einsum('kl,kl->kl',B,B), T) +a.T @ np.einsum('ij,ij->ij',A,A) @T -(2 * np.einsum('ab,ab->b', A @T @B, T)))

    


    fused_lgd1 = (1-alpha) * np.einsum('ik,ik->i', M,T)
    fused_lgd2 = (1-alpha) * np.einsum('ik,ik->k', M,T)

    lgd1 = geo_lgd1 + fused_lgd1
    lgd2 = geo_lgd2 + fused_lgd2

    return lgd1, lgd2




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
        pre_lgd = np.sum(ipdm_pp_w, axis = 0) #unclear which axis this should be
        eccentricity = pre_lgd**(1/p)

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
            




















