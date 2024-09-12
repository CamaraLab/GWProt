"""
This file contains data that can be used for FGW in addition to the usual isoelectric point

There are two ways to do it
- each residue is assigned a value and the fused component is the difference in values
    -pI
    -solvent accessible surface area
    -hydrophobicity
- there is a precomputed metric of distances
    -BLOSUM matrices
    -Grantham matrix

These are nested dicts of the form

Data[aa1][aa2] = d(aa1,aa2)
where aa1, aa2 are the one-letter amino acid codes

all are modified so as to be actual metrics when raw == False

raw == True  for those without corrections

the pI ignores the termini 

"""

import  re
import os
import math
import numpy as np
 
import Bio.PDB
from Bio import PDB, SeqIO
import blosum
import scipy.sparse.csgraph





canonical_aas = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

three_to_one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

missing_aa = '*'
def rescale_dictionary(d, rescaler):
    new_dict = {}
    for k in canonical_aas:
        new_dict[k] = {l:rescaler(d[l][k]) for l in canonical_aas  }
    return new_dict

def add_extra_aa(d,default = 'max'):
    # this one changes the dict to only have the 20 canonical aas and '*'
    # the distance to missing_aa is d, either the maximum or mean of the distances among the 20
    M = np.zeros((21,21))
    
    for i in range(20):
        for j in range(i):
            M[i,j] = d[canonical_aas[i]][canonical_aas[j]]
            M[j,i] = d[canonical_aas[i]][canonical_aas[j]]

    if default == 'max':
        de = np.max(M)
    if default =='avg':
        de = np.mean(M)
    if default == 'med':
        de = np.median(M)
    
    for i in range(20):
        M[i,20] = de
        M[20,i] = de

    new_dict = {}
    for i in range(20):
        i_dict = {}
        for j in range(20):
            i_dict[canonical_aas[j]] = M[i,j]
        i_dict[missing_aa] = de
        new_dict[canonical_aas[i]] = i_dict

    x_dict = { a: de for a in canonical_aas}

    new_dict[missing_aa] = x_dict
    new_dict[missing_aa][missing_aa] = de
    return new_dict
        

def correct_dictionary(d):
    # takes in a "sorta" - distance matrix in distance form, which may not satisfy the triangle inequality and returns one which does
    #assumes it's symmetric, undefined result if not
    keys = list(d.keys())
    N = len(keys)
    for k in keys:
        assert keys == list(d[k].keys())
    M = np.zeros((N,N))
    for i in range(N):
        for j in range(i):
            M[i,j] = d[keys[i]][keys[j]]
            M[j,i] = d[keys[i]][keys[j]]
    M = np.maximum( np.zeros((N,N)) , M)
    
    MM = scipy.sparse.csgraph.floyd_warshall(M, directed = False)
    new_dict = {}
    for i in range(N):
        i_dict = {}
        for j in range(N):
            i_dict[keys[j]] = MM[i,j]
        new_dict[keys[i]] = i_dict
    return new_dict

def dict_to_mat(d):
    keys = list(d.keys())
    N = len(keys)
    for k in keys:
        assert keys == list(d[k].keys())
    M = np.zeros((N,N))
    for i in range(N):
        for j in range(i):
            M[i,j] = d[keys[i]][keys[j]]
            M[j,i] = d[keys[i]][keys[j]]
    return M

def get_BLOSUM(n = 62, lamd = 1 , raw = False, default = 'med'):
    # n = 45,50,62,80,90
    B = blosum.BLOSUM(n)
    BB_raw = add_extra_aa(rescale_dictionary(B, lambda x : math.exp( -1 * lamd * x)), default = default)

    if raw:
        return BB_raw
    else:
        return correct_dictionary(BB_raw)

def get_Grantham(raw = False, default = 'med'):
    

    # Grantham 
    # https://www.jstor.org/stable/1739007?seq=1
    # https://github.com/maialab/grantham/blob/master/data-raw/grantham_distance_matrix.csv
    
    G = np.array([[0,110,145,74,58,99,124,56,142,155,144,112,89,68,46,121,65,80,135,177],
    [110,0,102,103,71,112,96,125,97,97,77,180,29,43,86,26,96,54,91,101],
    [145,102,0,98,92,96,32,138,5,22,36,198,99,113,153,107,172,138,15,61],
    [74,103,98,0,38,27,68,42,95,114,110,169,77,76,91,103,108,93,87,147],
    [58,71,92,38,0,58,69,59,89,103,92,149,47,42,65,78,85,65,81,128],
    [99,112,96,27,58,0,64,60,94,113,112,195,86,91,111,106,126,107,84,148],
    [124,96,32,68,69,64,0,109,29,50,55,192,84,96,133,97,152,121,21,88],
    [56,125,138,42,59,60,109,0,135,153,147,159,98,87,80,127,94,98,127,184],
    [142,97,5,95,89,94,29,135,0,21,33,198,94,109,149,102,168,134,10,61],
    [155,97,22,114,103,113,50,153,21,0,22,205,100,116,158,102,177,140,28,40],
    [144,77,36,110,92,112,55,147,33,22,0,194,83,99,143,85,160,122,36,37],
    [112,180,198,169,149,195,192,159,198,205,194,0,174,154,139,202,154,170,196,215],
    [89,29,99,77,47,86,84,98,94,100,83,174,0,24,68,32,81,40,87,115],
    [68,43,113,76,42,91,96,87,109,116,99,154,24,0,46,53,61,29,101,130],
    [46,86,153,91,65,111,133,80,149,158,143,139,68,46,0,94,23,42,142,174],
    [121,26,107,103,78,106,97,127,102,102,85,202,32,53,94,0,101,56,95,110],
    [65,96,172,108,85,126,152,94,168,177,160,154,81,61,23,101,0,45,160,181],
    [80,54,138,93,65,107,121,98,134,140,122,170,40,29,42,56,45,0,126,152],
    [135,91,15,87,81,84,21,127,10,28,36,196,87,101,142,95,160,126,0,67],
    [177,101,61,147,128,148,88,184,61,40,37,215,115,130,174,110,181,152,67,0]])
    
    # indices for the Grantham matrix:
    aas = "Ser Arg Leu	Pro	Thr	Ala	Val	Gly	Ile	Phe	Tyr	Cys	His	Gln	Asn	Lys	Asp	Glu	Met	Trp".upper().split()
    
    Grantham_raw = {}
    for k in three_to_one.keys():
        k_dict = {}
    
        for l in three_to_one.keys():
            i = aas.index(k)
            j = aas.index(l)
            k_dict[three_to_one[l]] = G[i,j]
        
        Grantham_raw[three_to_one[k]] = k_dict
    
    Grantham_raw = add_extra_aa(Grantham_raw, default = default)
    if raw:
        return Grantham_raw
    else:
        return correct_dictionary(Grantham_raw)


 

def get_pI():

    solomon = {"K" :  10.5,
    "R" :  12.5,
    "H" :  6,
    "D" :  3.9,
    "E" :  4.3,
    "C" :  8.3,
    "Y" :  10.1}
    
    for c in canonical_aas + [missing_aa]:
        if c not in solomon.keys():
            solomon[c] = 7
    
    pI_dict = {}

    for k in solomon.keys():
        k_dict = {k2: abs(solomon[k] -solomon[k2]  ) for k2 in solomon.keys()  }
        pI_dict[k] = k_dict
    return pI_dict
    #ignored:
    #C_term_pK = 2.4
    #N_term_pK =  9.6
    

def get_hydrophobicity():
    # Source: http://us.expasy.org/tools/pscale/Hphob.Eisenberg.html
# Amino acid scale: Normalized consensus hydrophobicity scale
# Author(s): Eisenberg D., Schwarz E., Komarony M., Wall R.
# Reference: J. Mol. Biol. 179:125-142 (1984)
    HH = {'C': 0.29,
 'D': -0.9,
 'S': -0.18,
 'Q': -0.85,
 'K': -1.5,
 'I': 1.38,
 'P': 0.12,
 'T': -0.05,
 'F': 1.19,
 'N': -0.78,
 'G': 0.48,
 'H': -0.4,
 'L': 1.06,
 'R': -2.53,
 'W': 0.81,
 'A': 0.62,
 'V': 1.08,
 'E': -0.74,
 'Y': 0.26,
 'M': 0.64,
 missing_aa: 0}
    HH_dict = {}

    for k in HH.keys():
        k_dict = {k2: abs(HH[k] -HH[k2]  ) for k2 in HH.keys()  }
        HH_dict[k] = k_dict
    return HH_dict


# deleted bc not compatible with pymol 3
# def get_SASA(file):
#     #returns solvent accessible surface area according to pymol
#     # this seems to work, but could be annoying to integrate
    

#     from pymol import cmd
#     import pymol
    
    
#     pymol.finish_launching(['pymol', '-pc'])
#     cmd.delete('all')
#     cmd.load(file, 'prot1')
#     residue_numbers = list(set(atom.resi_number for atom in cmd.get_model('prot1').atom))
    
    
#     SAS = [cmd.get_area(f'prot1 and resi {i}') for i in residue_numbers ]
#     return SAS


