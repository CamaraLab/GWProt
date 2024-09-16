import os
import numpy as np
import Bio.PDB
import Bio.SVDSuperimposer
import warnings
import ot
import math

from numpy import dot, transpose, sqrt
from numpy.linalg import svd, det

"""
https://biopython.org/docs/1.75/api/Bio.PDB.Superimposer.html#Bio.PDB.Superimposer.Superimposer

https://biopython.org/docs/1.75/api/Bio.SVDSuperimposer.html

source:
https://github.com/biopython/biopython/blob/master/Bio/SVDSuperimposer/__init__.py

"""

def Bio_RMSD(X,Y):
    #wrapper for the Bio.SVDSuperimposer
    sup = Bio.SVDSuperimposer.SVDSuperimposer()
    sup.set(X,Y)
    sup.run()
    rot, tran = sup.get_rotran()
    #print(rot,tran)
    return rot, tran

# def transform(X, rot = np.identity(3), trans = np.zeros((1,3))):
#     # applies the rotation matrix to X, then translates
#     assert rot.shape == (3,3)
#     assert trans.shape ==(3,)
#     assert X.shape[1] ==3
#     #( (rot @ X.T).T + trans ==  dot(X, rot) + trans).all()
#     return  (rot @ X.T).T + trans

#THIS NEEDS more TESTING
#no currently known bugs though

def weighted_RMSD(X:np.array ,Y:np.array, T:np.array):
    """
    This method uses the Kabsch algorithm to find a rigid, orientation-preserving transformation that minimizes weighted RMSD.

    Explicitly it finds a special orthogonal matrix S which minimizes

    sum_{i,j} |(x_i - x')- S(y_j - y')|^2 * T_{i,j} 

    where x' is the weighted mean of the x_i and y' is the weighted mean of the y'.

    Note - in general there may not be a unique solution matrix S which minimizes this. 

    :pararm X: A np.array with of (n,3) representing the n row vectors in R^3 defining X
    :pararm Y: A np.array with of (m,3) representing the m row vectors in R^3 defining Y
    :param T: A np.array of shape (n,m) representing the weights for the alignment. Its entries must be non-negative.
    So T[i,j] defines how strongly the distance between X[i,:] and Y[j,:] should be weighted in the minimization.
    :return: -y_mean, rot, x_mean ; where rot is a 3x3 matrix.


    Thus (Y - y_mean) @ rot + x_mean and X should be superimposed.
    Note that if n == m and T is the identity matrix, this is just the usual Kabsch algorithm for minimizing RMSD between X and Y.

    """
    
    #return pretranslation, rotation, posttranslation
    # assumed weights sum to 1



    # Y is the mobile one here, so
    # X ~ (Y+ pre) @ rot +post
        
        #code adapted from the Bio.SVDSuperimposer package:
        
    # Copyright (C) 2002, Thomas Hamelryck (thamelry@vub.ac.be)
    #
    # This file is part of the Biopython distribution and governed by your
    # choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
    # Please see the LICENSE file that should have been included as part of this
    # package.

    weights = T
    assert len(X.shape) == len(Y.shape) == len(weights.shape)

    assert X.shape[0] == weights.shape[0]
    assert Y.shape[0] == weights.shape[1]
    assert X.shape[1] ==3
    assert Y.shape[1] ==3
    assert (weights >=0 ).all()

    
    new_X = []
    new_Y = []
    new_weights = []
    counter = 0
    for i in range(weights.shape[0]):
        for j in range(weights.shape[1]):
            if weights[i,j] == 0:
                continue
            else:
                new_X.append(X[i,:])
                new_Y.append(Y[j,:])
                new_weights.append(weights[i,j])
                counter +=1
        
    new_X = np.stack(new_X)
    new_Y = np.stack(new_Y)
    W = np.diag(transpose(new_weights))#diagonal weight matrix
    # print('W.shape',W.shape)
    new_weights = np.stack(new_weights)[:, np.newaxis] #??? seems fishy

    

    #compute weighted means and translations
    # x_mean = np.sum(new_X * np.resize(new_weights, (counter,3)),axis = 0)  
    # y_mean = np.sum(new_Y * np.resize(new_weights, (counter,3)),axis = 0)  
    x_mean = np.sum(X * np.resize(np.sum(weights, axis = 1) , (X.shape[0],3)),axis = 0)/np.sum(weights)
    y_mean = np.sum(Y * np.resize(np.sum(weights, axis = 0) , (Y.shape[0],3)),axis = 0)/np.sum(weights)
#y_mean = np.sum(Y * np.resize(np.sum(weights, axis = 1) , (Y.shape[0],3)),axis = 0)  
    


    new_X = new_X - np.resize(x_mean, (counter,3))
    new_Y = new_Y - np.resize(y_mean, (counter,3))


    A = np.transpose(new_X)@ W @ new_Y
    u, d, vt = np.linalg.svd(A)
    rot = np.dot(np.transpose(vt), np.transpose(u))
    # check if we have found a reflection
    if det(rot) < 0:
        vt[2] = -vt[2]
        rot = dot(transpose(vt), transpose(u)) # is this supposed to be transposed ???
    #tran = x_mean - dot(y_mean, rot)
    return  -1 * y_mean, rot, x_mean
    # 



def get_raw_RMSD_stresses(X,Y):
    """
    This takes two point clouds of the same shape and computes the RMSD as well as the stresses
    assumes they are already superimposed
    """

    assert X.shape == Y.shape
    assert len(X.shape) == 2
    diff = X-Y
    sq_dists = np.sum(diff**2, axis = 1)
    mean_sq_dists = sq_dists / np.size(sq_dists)
    RMSD = math.sqrt(mean_sq_dists)
    return RMSD, mean_sq_dists

def get_RMSD_stresses(X,Y):
    """
    first applies the transformation then does get_raw_RMSD_stresses
    """
    pre, rot, post = weighted_RMSD(X,Y, np.identity(X.shape[0]))
    Y2 = (Y+ pre) @ rot +post
    return get_raw_RMSD_stresses(X,Y2)








def _pymol_transform( pretrans, rot, posttrans, Bio_format = True):
    # helper that changes this into the matrix format specific to pymol

    #Bio_format: row vectors, right multiplication. This is how the RMSD alignments do it
    # not Bio_format: column vectors, left multiplication
    warnings.warn('deprecation - this may use different target/mobile conventions than other functions')
    if not Bio_format:
        ll = list(rot[0, 0:3]) + [posttrans[0]] + list(rot[1, 0:3]) + [posttrans[1]] + list(rot[2, 0:3]) + [posttrans[2]] + [pretrans[0],pretrans[1],pretrans[2],1]
    if Bio_format:
        return(_pymol_transform(rot = rot.T, pretrans = pretrans, posttrans = posttrans, Bio_format = False))
    return(ll)


    





