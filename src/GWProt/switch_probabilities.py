import math
import numpy as np
import numpy.typing as npt
import itertools as it
from scipy.spatial.distance import *
from scipy.signal import *



import warnings



import sparse

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D
import matplotlib.patches as patches 

import mpl_toolkits.axisartist.floating_axes as floating_axes


def visualize_switch_probabibilities(A: np.array) -> None:
    """

    This method displays a switch probability matrix with matplotlib.

    :param A: The switch probability matrix to display; only the upper triangular part of it is used.

    """
    N = A.shape[0]    
    fig = plt.figure()
    plot_extents = N, 0, 0, N
    helper = floating_axes.GridHelperCurveLinear(Affine2D().rotate_deg(-45) + Affine2D().translate(tx = -math.sqrt(2)/2, ty = 0), plot_extents)
    ax = floating_axes.FloatingSubplot(fig, 111, grid_helper=helper)
    fig.add_subplot(ax)
    mask = np.triu(np.ones(A.shape), k=0)
    masked_arr = np.where(mask, A, np.nan)
    ax.imshow(masked_arr, cmap='hot', interpolation='none', transform=Affine2D(matrix =np.array([[0.70710678, 0.70710678, -0.        ],
           [ 0.70710678, -0.70710678,  0.        ],
           [ 0.        ,  0.        ,  1.        ]])) + ax.transData)
    
    ax.axis["bottom"].set_visible(False)
    ax.axis["left"].set_visible(False)
    ax.axis["top"].set_visible(True)
    ax.axis["right"].set_visible(True)
    
    ax.axis["top"].major_ticks.set_tick_out(True)
    ax.axis["right"].major_ticks.set_tick_out(True)
    
    
    
    
    
    ax.text(s = 'First Residue', x = N/(2*math.sqrt(2)) - 20, y = N/(2*math.sqrt(2)), rotation = 45, ha = 'center')
    ax.text(s = 'Second Residue', x = N*math.sqrt(2) - N/(2*math.sqrt(2)) + 20 , y = N/(2*math.sqrt(2)), rotation = -45, ha = 'center')
    
    ax.text(s = 'Residue Indices', x = N/math.sqrt(2)  , y = -N/4, ha = 'center')
    
    
    
    
    ax.axis["right"].major_ticks.set_visible(False)
    ax.axis["top"].major_ticks.set_visible(False)
    
    rect = patches.Rectangle((-2, 0), 2*N + 2, -1 * N, linewidth=1, edgecolor='black', facecolor="w" ) 
    ax.add_patch(rect)     
    
    for  i in range(N//100 + 1):
        ax.text(s = str(i*100), x = i*100*math.sqrt(2), y = -30, ha ='center', va = 'top')
    
    for i in range(N//10 +1):
        ax.text(s = '|', x = i*10*math.sqrt(2), y = 0, ha = 'center', va= 'top', weight = 'extra bold', fontsize = 1000/N+3)
        if i %10 ==0:
            ax.text(s = '|', x = i*10*math.sqrt(2), y = 0, ha = 'center', va= 'top', weight = 'extra bold', fontsize = 2000/N+3)
    
    plt.show()

def get_switch_probabilities(T: np.array, prot_num: int = 0) -> np.array:
    """
    Calculates the probability that the order of two residues are switched or not when the transport plan is applied.
    This can be used to detect circular permutations between two proteins. 

    :param T: The transport plan to use
    :param prot_num: Which protein to use, 0 uses the 0th axis of ``T``, 1 uses the 1st axis.
    :return: A square np.array whose *ij*th entry is the probability that residues *i* and *j* are kept in the same order. 
    
    """
    
    if prot_num == 1:
        return get_switch_prob(T.T, prot_num = 0)

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




def preprocess(A: np.array) -> np.array:
    """

    Processes a switch probability matrix for applying ``max_rectangle_diagonal``

    :param A: A switch probability matrix
    :return: A processed switch probability matrix

    """
    kernel = np.ones((5,5))
    mat2 = oaconvolve(A, kernel, mode = 'same')/np.sum(kernel)
    return np.triu((1 - mat2 + 0.6)).astype(int)



def _histogram_area_helper_diag(histogram,j, min = 0):
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
            if area >= max_area and ((index -1  +histogram[top_of_stack] -1 ==j) or (((stack[-1] +1) if stack else 0) -histogram[top_of_stack] +1 ==j)) and histogram[top_of_stack] >= min  and ((index - stack[-1] - 1) if stack else index) >= min:
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
        if area >= max_area and ((index -1  +histogram[top_of_stack] -1 ==j) or (((stack[-1] +1) if stack else 0) -histogram[top_of_stack] +1 ==j)) and histogram[top_of_stack] >= min  and ((index - stack[-1] - 1) if stack else index) >= min:
            max_area = area
            #print(j, area, index, (stack[-1] +1) if stack else 0, top_of_stack) #debugging
            left, right  =   (stack[-1] +1) if stack else 0,index -1
            height = histogram[top_of_stack]
            
            #print('T2,',  left, right,index -1 , area, stack)
    # Return maximum area under the given histogram 
    return max_area ,left , right, height


def max_rectangle_diagonal(A : np.array, min: int = 0) -> tuple[int, tuple[int,int,int,int]]: 
    """

    This method finds the area and coordinates of the largest rectange in an array containing only nonzero entries, with a corner on the main diagonal.
    
    :param A: The array
    :param min: Rectangles whose width or height are below ``min`` are not considered
    :return: This returns the maximal area and the coordinates of the rectangle in the tuple ``(max_area, (left,  right, bottom, top ))`` 
    
    """
    
    mat = (A).astype(int)
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
        area , bottom, top, height = _histogram_area_helper_diag(histogram = l,j = j, min = min)
        # if area > best_area and (top +height -1 ==j) or (bottom -height +1 ==j):
        if area > best_area and (top-bottom + 1) >= min and height >= min:
            #print(j, bottom, top, height)


            best_area = area
            if j >= bottom: #top half
                best_coords = (  int(j - height+1), int(j), int(bottom), int(top) ) 
                #print('top', best_coords)
            elif top >= j: #bottom half
                best_coords = ( int(j) , int(j + height-1), int(bottom), int(top) ) 
                #print('bot', best_coords)
    return best_area , best_coords




    