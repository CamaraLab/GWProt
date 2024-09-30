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
from scipy.signal import *
import multiprocessing
import multiprocess

import pandas as pd
import cajal

import Bio.PDB
from typing import Iterator, Iterable, Optional, TypeVar, Generic



import warnings

import sys
sys.path.insert(0,'../PGC020.a12/src')

import GW_scripts
import read_pdb
import FGW_protein
import IdInit
import sparse

    @staticmethod
    def get_switch_prob_sparse(T: np.array, prot_num: int = 0):
        """
        Calculates the probability that the order of two residues are switched or not when the transport plan is applied.
        This can be used to detect circular permutations between two proteins. This uses sparse matrices so may cause problems if T has many non-zer entries.
        :param T: The transport plan to use
        :param prot_num: Which protein to use, 0 uses the 0th axis of 'T', 1 uses the 1st axis.
        :return: A square np.array whose ijth entry is the probability that residues i and j are kept in the same order. 
        """
        
        if prot_num == 1:
            return FGW_protein.get_switch_prob_sparse(T.T, prot_num = 0)

        if np.count_nonzero(np.sum(T, axis = 1) ==0) > 0:
            raise ValueError('T has a zero row or column')
        T_mod = T/ (np.sum(T, axis = 1)[np.newaxis]).T
        TT_mod = T_mod.T 
        if np.count_nonzero(T) >= 5000:
            warnings.warn('input has over 5,000 nonzero entries, may use too much RAM and crash')
            
        
        T_sparse = sparse.COO.from_numpy(T_mod)
        TT_sparse = sparse.COO.from_numpy(TT_mod)
        #print(T_sparse.shape)
        #print(TT_sparse.Ashape)
        
       # sparse_big_one= sparse.einsum( 'il,jk-> klij'  ,T_sparse, TT_sparse) #this one is probably wrong
        sparse_big_one= sparse.einsum( 'il,kj-> ijkl'  ,T_sparse, TT_sparse)
        
        #  goal:
        #  (i,j,k,l)th entry is the ijth entry of the matrix given by the kth colum and lth column transposed
        L = sparse.tril(sparse_big_one, -1)
        P1 = sparse.COO.todense(sparse.COO.sum(L, (2,3)))
        
        # U = sparse.triu(sparse_big_one, 1)
        # P2 = sparse.COO.todense(sparse.COO.sum(U, (2,3))) #this is just the transpose of P1, up to floating point inaccuracies

        return P1



def messiness_score(mat, normalize = True):
    #very basic convolutional filter for edge detection
    # oaconvolve not yet tested
    
    e = oaconvolve(mat, np.array([[0,-1,0],[-1,4,-1],[0,-1,0]]), mode = 'same')
    edge_score = np.sum(np.abs(e))
    
    # edge_mat = 4*mat[1:-1, 1: -1]
    
    # edge_mat -= mat[1:-1, 0:-2]
    # edge_mat -= mat[1:-1, 2:]
    # edge_mat -= mat[0:-2,1:-1]
    # edge_mat -= mat[2:,1:-1]
    # edge_mat = np.abs(edge_mat)
    # edge_score = np.sum(edge_mat)
    
    # now consider values not within epsilon of 0 or 1
    # max = np.max(mat)
    # m = lambda x : (0.05*max < x < 0.95 *max)
    # color_score = np.count_nonzero(np.vectorize(m)(mat) )   # doesn't seem as useful
    
    if normalize:
        return edge_score/mat.size 
    else:
        return edge_score

def preprocess_smooth(mat, kernel = np.ones((5,5))):
#not sure what other kernel to use
	return oaconvolve(mat, kernel, mode = 'same')/np.sum(kernel)

def histogram_area_helper(histogram, min = 0):
    #takes in an int list and returns the largest area rectangle it contains an the start and end indices of it
    # based on:
    #https://www.geeksforgeeks ,org/largest-rectangular-area-in-a-histogram-using-stack/
    # This code is contributed 
# by Jinay Shah 

    # Python3 program to find maximum 
# rectangular area in linear time 



    # This function calculates maximum 
    # rectangular area under given 
    # histogram with n bars 

    # Create an empty stack. The stack 
    # holds indexes of histogram[] list. 
    # The bars stored in the stack are 
    # always in increasing order of 
    # their heights. 
    stack = list() 

    max_area = 0 # Initialize max area 
    left = -1
    right = -1
    height = 0

    # Run through all bars of given histogram 
    index = 0
    while index < len(histogram): 
        # If this bar is higher than the bar on top  stack, push it to stack 
        if (not stack) or (histogram[stack[-1]] <= histogram[index]): 
            stack.append(index) 
            index += 1
        # If this bar is lower than top of stack, then calculate area of rectangle with 
        # stack top as the smallest (or minimum  height) bar.'i' is 'right index' for 
        # the top and element before top in stack  is 'left index' 
        else: 
            # pop the top 
            top_of_stack = stack.pop() 
            # Calculate the area with histogram[top_of_stack] stack as smallest bar 
            area = (histogram[top_of_stack] *
                    ((index - stack[-1] - 1) if stack else index)) 

            #print('T1,',   (stack[-1] +1) if stack else 0,index -1, area, stack)
            
            # update max area, if needed 
            if area >= max_area and histogram[top_of_stack] > min  and (index - stack[-1] - 1 if stack else index) > min:
                max_area = area
                left, right  =   (stack[-1] +1) if stack else 0,index -1
                height = histogram[top_of_stack]
            # elif area >= max_area:
            #     print(histogram[top_of_stack] , (index - stack[-1] - 1 if stack else index))

    # Now pop the remaining bars from stack and calculate area with every popped bar as the smallest bar 
    while stack: 
        #print(left,right)
        # pop the top 
        top_of_stack = stack.pop() 
        # Calculate the area with histogram[top_of_stack] stack as smallest bar 
        area = (histogram[top_of_stack] *
                ((index - stack[-1] - 1) 
                if stack else index)) 
        # update max area, if needed 
        if area >= max_area and histogram[top_of_stack] > min and ((index - stack[-1] - 1) if stack else index) > min:
            max_area = area
            left, right  =   (stack[-1] +1) if stack else 0,index -1
            height = histogram[top_of_stack]
            
            #print('T2,',  left, right,index -1 , area, stack)
        # elif area >= max_area:
        #     print(histogram[top_of_stack] , (index - stack[-1] - 1 if stack else index))
    # Return maximum area under the given histogram 
    return max_area ,left , right, height


def max_rectangle(Mat, min=0): 
    mat = (Mat+0.5).astype(int)
    #https://drdobbs.com/database/the-maximal-rectangle-problem/184410529
    m, n = mat.shape

    cache_mat = np.zeros(mat.shape)
    
    for i in range(m):
        if mat[i,-1] :
            cache_mat[i,-1] = 1
        for j in range(n-2,-1,-1):
            if mat[i,j]:
                cache_mat[i,j] = cache_mat[i,j+1] +1

    #print(cache_mat)

    best_area = 0
    best_coords = (0,n,0,m) # left, right, bottom, top , inclusive
    for j in range(n):
        
        l = list(cache_mat[:,j].T)
        area , bottom, top, height = histogram_area_helper(l, min = min)
        #print(j,area,  height, (top-bottom+1))
        if area > best_area  and (top-bottom + 1) > min and height > min:
            best_area = area
            best_coords = ( j , j + height-1, bottom, top ) 
    return best_area , best_coords

def histogram_area_helper_diag(histogram,j, min = 0):
    #takes in an int list and returns the largest area rectangle it touching the main diagonal contains an the start and end indices of it
    #j is the column index

    stack = list() 

    max_area = 0 # Initialize max area 
    left = None
    right = None
    height = 0

    # Run through all bars of given histogram 
    index = 0
    while index < len(histogram): 
        # If this bar is higher than the bar on top  stack, push it to stack 
        if (not stack) or (histogram[stack[-1]] <= histogram[index]): 
            stack.append(index) 
            index += 1
        # If this bar is lower than top of stack, then calculate area of rectangle with 
        # stack top as the smallest (or minimum  height) bar.'i' is 'right index' for 
        # the top and element before top in stack  is 'left index' 
        else: 
            # pop the top 
            top_of_stack = stack.pop() 
            # Calculate the area with histogram[top_of_stack] stack as smallest bar 
            area = (histogram[top_of_stack] *
                    ((index - stack[-1] - 1) if stack else index)) 

            #print('T1,',   (stack[-1] +1) if stack else 0,index -1, area, stack)
            
            # update max area, if needed 
            if area >= max_area and ((index -1  +histogram[top_of_stack] -1 ==j) or (((stack[-1] +1) if stack else 0) -histogram[top_of_stack] +1 ==j)) and histogram[top_of_stack] > min  and ((index - stack[-1] - 1) if stack else index) > min:
                max_area = area
                #print(j, area, index, (stack[-1] +1) if stack else 0, top_of_stack) #debugging
                left, right  =   (stack[-1] +1) if stack else 0,index -1
                height = histogram[top_of_stack]

    # Now pop the remaining bars from stack and calculate area with every popped bar as the smallest bar 
    while stack: 
        #print(left,right)
        # pop the top 
        top_of_stack = stack.pop() 
        # Calculate the area with histogram[top_of_stack] stack as smallest bar 
        area = (histogram[top_of_stack] *
                ((index - stack[-1] - 1) 
                if stack else index)) 
        # update max area, if needed 
        if area >= max_area and ((index -1  +histogram[top_of_stack] -1 ==j) or (((stack[-1] +1) if stack else 0) -histogram[top_of_stack] +1 ==j)) and histogram[top_of_stack] > min  and ((index - stack[-1] - 1) if stack else index) > min:
            max_area = area
            #print(j, area, index, (stack[-1] +1) if stack else 0, top_of_stack) #debugging
            left, right  =   (stack[-1] +1) if stack else 0,index -1
            height = histogram[top_of_stack]
            
            #print('T2,',  left, right,index -1 , area, stack)
    # Return maximum area under the given histogram 
    return max_area ,left , right, height


def max_rectangle_diagonal(Mat, min = 0): 
#min = minimum width or height of rectangle to allow
    
    mat = (Mat + 0.5).astype(int)
    #https://drdobbs.com/database/the-maximal-rectangle-problem/184410529
    assert mat.shape[0] == mat.shape[1]
    n = mat.shape[0]

    #create and fill cache
    cache_mat = np.zeros(mat.shape)

    #first initialize diagonal
    for i in range(n):
        cache_mat[i,i] = int(mat[i,i])

    #top half
    for i in range(n):
        for j in range(i,n):
            if mat[i,j]:
                cache_mat[i,j] = cache_mat[i,j-1] +1 
    #bottom half
    for i in range(1,n):
        for j in range(i-1,-1,-1):
            if mat[i,j]:
                cache_mat[i,j] = cache_mat[i,j+1] +1 
    best_area = 0
    best_coords = (-1,-1,-1,-1) # left, right, bottom, top , inclusive
    for j in range(n):
        
        l = list(cache_mat[:,j].T)
        
        #print(l)
        area , bottom, top, height = histogram_area_helper_diag(histogram = l,j = j, min = min)
        # if area > best_area and (top +height -1 ==j) or (bottom -height +1 ==j):
        if area > best_area and (top-bottom + 1) > min and height > min:
            #print(j, bottom, top, height)


            best_area = area
            if j >= bottom: #top half
                best_coords = (  j - height+1, j, bottom, top ) 
                #print('top', best_coords)
            elif top >= j: #bottom half
                best_coords = ( j , j + height-1, bottom, top ) 
                #print('bot', best_coords)
    return best_area , best_coords




    