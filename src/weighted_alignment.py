import os
import numpy as np
import Bio.PDB
import Bio.SVDSuperimposer
import warnings
import ot

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

def transform(X, rot = np.identity(3), trans = np.zeros((1,3))):
    # applies the rotation matrix to X, then translates
    assert rot.shape == (3,3)
    assert trans.shape ==(3,)
    assert X.shape[1] ==3
    #( (rot @ X.T).T + trans ==  dot(X, rot) + trans).all()
    return  (rot @ X.T).T + trans

#THIS NEEDS more TESTING
#no currently known bugs though

def weighted_RMSD(X,Y, weights):
    
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
    
    assert X.shape[0] == weights.shape[0]
    assert Y.shape[0] == weights.shape[1]
    assert X.shape[1] ==3
    assert Y.shape[1] ==3

    # print('X.shape',X.shape)
    # print('Y.shape',Y.shape)
    
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


"""
transform based on matrix and vector in pymol:
"""




def pymol_transform( pretrans, rot, posttrans, Bio_format = True):
    #Bio_format: row vectors, right multiplication. This is how the RMSD alignments do it
    # not Bio_format: column vectors, left multiplication
    warnings.warn('deprecation - this may use different target/mobile conventions than other functions')
    if not Bio_format:
        ll = list(rot[0, 0:3]) + [posttrans[0]] + list(rot[1, 0:3]) + [posttrans[1]] + list(rot[2, 0:3]) + [posttrans[2]] + [pretrans[0],pretrans[1],pretrans[2],1]
    if Bio_format:
        return(pymol_transform(rot = rot.T, pretrans = pretrans, posttrans = posttrans, Bio_format = False))
    return(ll)


# 
def reform_transport_plan(T, coords1, coords2):
    # takes in two sets of coordinates with a transport plan
    # applies the rigid translation minimizing weighted RMSD according to the transport plan (of weights)
    # 
    n ,m = T.shape
    assert coords1.shape == (n,3)
    assert coords2.shape == (m,3)

    pretrans, rot, posttrans = weighted_RMSD( coords1 , coords2, T)
    shifted_coords2 = (coords2 + pretrans) @ rot + posttrans # this is the way
    D = ot.dist(coords1, shifted_coords2)
    #print(D.shape)
    a = np.ones(coords1.shape[0])/coords1.shape[0]
    b = np.ones(shifted_coords2.shape[0])/coords2.shape[0]

    TT = ot.emd(a,b, D)

    stress = np.einsum('ij,ij ->ij', D,T)
    cost = np.sum(stress)
    stress1 = np.sum(stress, axis = 1)
    stress2 = np.sum(stress, axis = 0)

    return TT , cost, stress1, stress2 #pretrans, rot, posttrans, shifted_coords2

    
    
    





