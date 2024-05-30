"""
read_pdb.py

A function that takes in a pdb file and makes it parseable for CAJAL
4/1/2024
"""
import numpy as np
import Bio.PDB.PDBParser
import numpy.typing as npt
import re
import statistics

def get_pdb_coords(
    # this takes in a pdb file and gets the coordinates of the CA atoms, outputting them as an array
    # code mostly by Patrick and Pablo
    filepath : str,
    n : int = np.inf #the max number of samples we take
    ): #returns array of shape (n,3) of 3D coords of the CA's chosen amino acids
    """
    filepath is a full filepath.
    Returns an (n, 3) ndarray of points, n <= N, each row is a point in Euclidean space.
    """
    
    pattern =  r'\/?([^\/]+)\.pdb$'

    cellname = re.findall(pattern,filepath )[0]
    parser = Bio.PDB.PDBParser() #aaaaahahhaha
    #print(cellname,filepath) #for debugging
    
    cell = parser.get_structure(cellname,filepath)
    #print(cell)
    model = next(cell.get_models())
    #print("model is")
    #print(model)
    residues = model.get_residues()
    #print(residues)
    coords = []
    for res in residues:

        het_flag = res.get_id()[0]
        #print(res.get_resname())
        if het_flag !=' ':
            continue
        #sort(res)
        atoms = res.get_atoms()  #N, CA, C ,  CB and O order varies
        N = next(atoms)
        assert N.get_id() == "N"
        N_coords = N.get_coord()
        N_vect = np.array(N_coords)
        
        CA = next(atoms)
        assert CA.get_id() == "CA"
        CA_coords = CA.get_coord()
        CA_vect = np.array(CA_coords)
        
        if True:
            coords.append(CA.get_coord())
            continue
        
        C = next(atoms)
        assert C.get_id() == "C"
        C_coords = C.get_coord()
        C_vect = np.array(C_coords)
        
        
        
        if res.get_resname() == "GLY": 
            #calculate the imaginary CB location
            #this has not yet been tested
            N_C_avg = (N_vect + C_vect) /2
        
            CBh = CA_vect - N_C_avg #component of CB in the plane of CA,C,N
        
        
            perp = np.cross( N_vect - CA_vect, C_vect - CA_vect) / np.linalg.norm(np.cross( N_vect - CA_vect, C_vect - CA_vect) )
        
            CBv =  perp* np.linalg.norm(C_vect - N_vect)/2
        
            CB_vect = CBv + CBh+ CA_coords
            CB_coords = list(CB_vect)
        
            #seems correct now
        else:
            CB = next(atoms)
            if CB.get_id() == "O" : 
                CB = next(atoms)
            assert CB.get_id() == "CB"
            CB_coords = CB.get_coord()
            CB_vect = np.array(CB_coords)   	
        
        if True:
            coords.append(CA.get_coord())
        
        del CA
        del res
        
        #del model
        #del cell
        
        
        
    arr=np.stack(coords)
    if len(coords)>n:
        new_indices=np.linspace(0, len(coords), num=n, endpoint=False, dtype=np.int_) # samples n of them, evenly spaced
        arr=arr[new_indices,:]
    #print(arr[:5]) #for debugging
    return arr #, residues, model, cell

def get_pdb_coords_pI(
    # this takes in a a protein pdb file and outputs 2 arrays
    # the first is a list of length <= n of coordinates, 
    # the second is of estimated isoelectric points of the oligopeptdes 
    filepath : str,
    n : int = np.inf, #the max number of samples we take
    median: bool = False
    ):
    
    location_type = "CA"
    pKSet = "solomon"
    pattern =  r'\/([^\/]+)\.pdb$' 
    
    cellname = re.findall(pattern,filepath )[0]
    parser = Bio.PDB.PDBParser() #aaaaahahhaha
    #print(cellname,filepath) #for debugging
    
    cell = parser.get_structure(cellname,filepath)
    #print(cell)
    model = next(cell.get_models())
    #print("model is")
    #print(model)
    residues = model.get_residues()
    #print(residues)
    coords = []
    res_list = ""
    
    for res in residues:
        if 'UNK' in res.get_resname():
            res_list += 'X' #unknown
        else:
            res_list += res_3_to_1_dict[res.get_resname() ] # adds the 
    
        het_flag = res.get_id()[0]
        #print(res.get_resname())
        if het_flag !=' ':
            continue
        #sort(res)
        atoms = res.get_atoms()  #N, CA, C ,  CB and O order varies
        N = next(atoms)
        assert N.get_id() == "N"
        N_coords = N.get_coord()
        N_vect = np.array(N_coords)
    
        CA = next(atoms)
        assert CA.get_id() == "CA"
        CA_coords = CA.get_coord()
        CA_vect = np.array(CA_coords)
    
        if location_type == "CA":
            coords.append(CA.get_coord())
            continue
    
        C = next(atoms)
        assert C.get_id() == "C"
        C_coords = C.get_coord()
        C_vect = np.array(C_coords)
    
        
    
        if res.get_resname() == "GLY": 
            #calculate the imaginary CB location
            #this has not yet been tested
            N_C_avg = (N_vect + C_vect) /2
    
            CBh = CA_vect - N_C_avg #component of CB in the plane of CA,C,N
    
    
            perp = np.cross( N_vect - CA_vect, C_vect - CA_vect) / np.linalg.norm(np.cross( N_vect - CA_vect, C_vect - CA_vect) )
    
            CBv =  perp* np.linalg.norm(C_vect - N_vect)/2
    
            CB_vect = CBv + CBh+ CA_coords
            CB_coords = list(CB_vect)
    
            #seems correct now
        else:
            CB = next(atoms)
            if CB.get_id() == "O" : 
                CB = next(atoms)
            assert CB.get_id() == "CB"
            CB_coords = CB.get_coord()
            CB_vect = np.array(CB_coords)   	
    
        if location_type == "CA":
            coords.append(CA.get_coord())
        elif location_type == "CB":
            coords.append(CB_coords)
        elif location_type == "virtual_center":
    
            Nproj = N_vect -  np.dot(N_vect - CA_vect , CB_vect - CA_vect )* (CB_vect - CA_vect)  #takes CA as origin
    
            direction= -1* Nproj/np.linalg.norm(Nproj)
    
            VC_vect = CA_vect + 2* np.linalg.norm(CB_vect - CA_vect) * direction
            coords.append(list(VC_vect))
            
        del CA
        del res
    
    #del model
    #del cell
    
    #print('test1')
    
    #arr=np.stack(coords)
    
    coord_segments = split_list(coords, n)
    res_segments = split_list(res_list,n)
    assert len(coord_segments) == len(res_segments)
    
    coord_list = []
    pI_segment_list = []
    #print('test2')
    #figure out how to deal with N, C terminii
    #assumes first is N, last is C
    
    
    for i in range(len(coord_segments)):
        e = np.stack(coord_segments[i]) 
    
        coord_list.append(np.mean(e,0) )
        
    if median:
        pI_segment_list.append(writeProtIepMedian(res_segments[0],  N_term =True))
        for seg in res_segments[1:-1]:
            pI_segment_list.append(writeProtIepMedian(seg))
    
        pI_segment_list.append(writeProtIepMedian(res_segments[-1], C_term =True))
    
    else:
        pI_segment_list.append(writeProtIep(res_segments[0],  N_term =True))
        for seg in res_segments[1:-1]:
            pI_segment_list.append(writeProtIep(seg))
    
        pI_segment_list.append(writeProtIep(res_segments[-1],  C_term =True))
    
    
    return coord_list, pI_segment_list

def split_list(l, n):
    #splits the list l into n (roughly) evenly sized smaller lists
    # if  len(l) % n != 0 then it tries to evenly space the different sized sublists
    # returns a list of lists
    
    length = len(l)
    if n <=1:
        return [l]
    if n >= length:
        return [[a] for a in l]
    
   
    return [l[int(i*length/n): int((i+1)*length/n) ]  for i in range(n) ]

"""
//    Sequence Manipulation Suite. A collection of simple JavaScript programs
//    for generating, formatting, and analyzing short DNA and protein
//    sequences.
//    Copyright (C) 2020 Paul Stothard stothard@ualberta.ca
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.
//

//Written by Paul Stothard, University of Alberta, Canada

https://github.com/paulstothard/sequence_manipulation_suite

adapted and changed to python by Elijah Gunther


Several other residue pKa value tables can be found here:
https://github.com/bigbio/pIR/blob/master/R/pkSets.R

"""


res_3_to_1_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def writeProtIepMedian(proteinSequence,  N_term: bool = False, C_term: bool = False):
    #instead uses the median of the charged/ionic residues + termini
    charges = []
    proteinSequence = str(proteinSequence).lower()


    N_term_pK =  9.6
    K_pK =  10.5
    R_pK =  12.5
    H_pK =  6
    D_pK =  3.9
    E_pK =  4.3
    C_pK =  8.3
    Y_pK =  10.1
    C_term_pK = 2.4

        
    if N_term:
        charges.append(N_term_pK)
    if C_term:
        charges.append(C_term_pK)
    for res in proteinSequence:
        if res == "k":
            charges.append(K_pK)
        if res == "r":
            charges.append(R_pK)
        if res == "h":
            charges.append(H_pK)    
        if res == "d":
            charges.append(D_pK)
        if res == "e":
            charges.append(E_pK)            
        if res == "c":
            charges.append(C_pK)
        if res == "y":
            charges.append(Y_pK)

    if len(charges) ==0:
        return 7
    else:
        return statistics.median(charges)




def writeProtIep(proteinSequence,  pKSet:str = "solomon", N_term: bool = False, C_term: bool = False):

    #THIS NEEDS TO BE CITED!!!
  #calculates pI of protein.
    # N_term, C_term are whether the N (resp C) terminus is in this part of the protein
    pH = 7.0
    step = 1 #originally 3.5
    charge = 0.0
    last_charge = 0.0
    
    proteinSequence = str(proteinSequence).lower()



    N_term_pK =  9.6
    K_pK =  10.5
    R_pK =  12.5
    H_pK =  6
    D_pK =  3.9
    E_pK =  4.3
    C_pK =  8.3
    Y_pK =  10.1
    C_term_pK = 2.4


    
    K_count = proteinSequence.count('k')
    R_count = proteinSequence.count('r')
    H_count = proteinSequence.count('h')
    D_count = proteinSequence.count('d')
    E_count = proteinSequence.count('e')
    C_count = proteinSequence.count('c')
    Y_count = proteinSequence.count('y')

  


    while True:
        charge = (K_count * partial_charge(K_pK, pH) +
        R_count * partial_charge(R_pK, pH) +
        H_count * partial_charge(H_pK, pH) -
        D_count * partial_charge(pH, D_pK) -
        E_count * partial_charge(pH, E_pK) -
        C_count * partial_charge(pH, C_pK) -
        Y_count * partial_charge(pH, Y_pK) )
        
        if N_term:
            charge += partial_charge(N_term_pK, pH)
        if C_term:
            charge -= partial_charge(pH, C_term_pK)
            
        #print(charge) #debugging
        #print(pH) #debugging
        
        
        if abs(charge) < 0.1:
        #round(charge,2) == round(last_charge*100,2):
            #print(charge, last_charge)
            break #check that this is same as in python
            
        # (charge.toFixed(2) == (last_charge * 100).toFixed(2))
        

        if charge > 0: 
            pH = pH + step
        else:
            pH = pH - step


        step = step*0.9 #originally 0.5

        last_charge = charge
  


    return pH


def partial_charge(first, second):
    charge = 10**( first - second)
    return charge / (charge + 1) 









