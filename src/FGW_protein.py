import os
os.environ['OPENBLAS_NUM_THREADS'] = '1' #should be first, before ot

import time
import re
import math
import numpy as np
import random
import ot
import statistics
import numpy.typing as npt
import itertools as it
from scipy.spatial.distance import *
from deprecated import deprecated

import Bio.PDB
from Bio import PDB, SeqIO
from typing import Callable


import IdInit
import GW_scripts
import read_pdb
import run_fasta36

from cajal import run_gw, qgw, gw_cython
"""
copied 4/1/2024
"""




def pH_median(l: list[float]) -> float:
    """
    algorithm for approximating the isoelectric point, takes the median after ignoring 7's
    returns 7 if the length is 0
    :param l: list of isoelectric points of residues

    :return: estimate of isoelectric point of polypeptide
    """


    
    ll = [p for p in l if p!=7]
    if len(ll) ==0:
        return 7
    else:
        return np.median(ll)

def pI_iter_alg(proteinSequence: str,
    N_term_count:int = 0,
    C_term_count:int = 0) -> float:
    """
    Estimates the isoelectric point of a polypeptide using the Henderson-Hasselbach equation
    The code is based on code in the Sequence Manipulation Suite, 
    and uses the pI values of residues and termini from Solomon’s Organic Chemistry, fifth edition

    Stothard P. The sequence manipulation suite: JavaScript programs for analyzing and formatting protein and DNA sequences. Biotechniques. 
    2000 Jun;28(6):1102, 1104. doi: 10.2144/00286ir01. PMID: 10868275.

    :param proteinSequence : string of the polypeptide's residues in one-letter symbols
    :param N_term_count : number of Nitrogen termini to include
    :param C_term_count : number of Carbon termini to include
    :return: the estimated isoelectric point
    """


        
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
    //    but WITHOUT ANY WARRANTY; without even the implied warranty of
    //    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    //    GNU General Public License for more details.
    //
    //    You should have received a copy of the GNU General Public License
    //    along with this program.  If not, see <https://www.gnu.org/licenses/>.
    //
    
    //Written by Paul Stothard, University of Alberta, Canada 
    rewritten in python and revised by Elijah Gunther (all errors are my own)
    
    https://github.com/paulstothard/sequence_manipulation_suite/blob/master/docs/scripts/protein_iep.js
    
    Stothard P. The sequence manipulation suite: JavaScript programs for analyzing and formatting protein and DNA sequences. Biotechniques. 
    2000 Jun;28(6):1102, 1104. doi: 10.2144/00286ir01. PMID: 10868275.



    Graham Solomons T.W., Fryhle C.B., Snyder S.A.. Solomons’ Organic Chemistry. 

    """
    #calculates pI of protein.
    # N_term, C_term are whether the N (resp C) terminus is in this part of the protein
    #and how many of them we include for the convolution

    pH = 7.0
    step = 1 #originally 3.5
    charge = 0.0
    last_charge = 0.0
    
    proteinSequence = str(proteinSequence).lower()
    
    #solomon
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
        Y_count * partial_charge(pH, Y_pK) +
        N_term_count * partial_charge(N_term_pK, pH) - 
        C_term_count * partial_charge(pH, C_term_pK)     )
    
    
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
    #helper method for pI_iter_alg
    charge = 10**( first - second)
    return charge / (charge + 1)



class FGW_protein:
    """
    This class contains everything needed to run GW and fused GW on proteins, as well as versions with distortion scaling and sequence alignment

    :param name: Simply for ease of use
    :param coords: The coordinates of the CA atoms of the protein, ordered sequentially
    :param fasta: A string in the format of a fasta file containing a header and the sequence of the protein. Sequence has length n
    :param ipdm: The intra-protein distance matrix of a protein. 
        The (i,j)th entry is the (possibly scaled) distance between residues i and j. This is mutable can can change if distortion scaling is used.
    :param scaled_flag: Records whether the ipdm is the exact distance between residues or if it has been scaled.

    """




    def __init__(self, name, fasta, pI_list,  coords = None, ipdm = None, scaled_flag = False ):
        #note - the fasta full text of the file, not the sequence
        #input validation
        assert not (coords is None and ipdm is None)
        if not coords is None:
            coords = np.array(coords)

            # print('type(coords)' ,type(coords)) #debugging
            # print('coords.shape' ,coords.shape)
            # print('len(pI_list)', len(pI_list))
            
            assert len(coords.shape) == 2
            assert coords.shape[1] == 3 
            assert coords.shape[0] == len(pI_list)
        
        if ipdm is not None:
            ipdm = np.array(ipdm)
            assert len(ipdm.shape) ==2
            assert ipdm.shape[0] == ipdm.shape[1]
            assert (ipdm == ipdm.T).all()

        fasta_header, fasta_seq = re.findall(string = fasta, pattern = r'^(>.+)\n([A-Z]*)$')[0]
        assert len(fasta_seq) ==len(pI_list)
        
        self.name = name
        self.fasta = fasta
        self.pI_list = pI_list
        self.coords = coords
        
        self.scaled_flag = scaled_flag #whether the ipdm has been scaled

        if ipdm is None:
            self.ipdm = squareform(pdist(self.coords))
        else:
            
            self.ipdm = ipdm
            
    def __eq__(self, other):
        """
        Compares the underlying fasta sequences (not the full fasta file), the pI_lists, the ipdms, and the coords if both are defined.
        This does NOT compare the names, scaled_flags, or fasta headers.
        """


        fasta_header1, fasta_seq1 = re.findall(string = self.fasta, pattern = r'^(>.+)\n([A-Z]*)$')[0]
        fasta_header2, fasta_seq2 = re.findall(string = other.fasta, pattern = r'^(>.+)\n([A-Z]*)$')[0]
        
        if self.coords is not None and other.coords is not None and (self.coords != other.coords).any():
            return False  
        return fasta_seq1 == fasta_seq2 and self.pI_list == other.pI_list and (self.ipdm == other.ipdm).all()
      

    def convolve_pIs_fasta(self, 
        kernel_list :list[int], 
        origin: int, 
        inplace :bool = False) -> list[float]:
        """
        This method applies a convolution process to the 'FGW_protein' object which smoothes out the isoelectic points associated
        to each residue by combining them with those of nearby residues. The intended use is that this could be applied before downsampling
        so that the isoelectric points of discarded residues is still preserved. That is done automatically with downsample_n(pI_combination = True),
        so this is most useful when applied before run_FGW_seq_aln() as that method discards unaligned residues.


        The convolution works as follows: for each residue we make a virtual oligopeptide of copies of that residue and its neighbors, 
        then use 'pI_iter_alg' to estimate the oligopeptide's isoelectric point. The number of copies is the entry of 'kernel_list', where the current residue is at position 'origin',
        The isoelectric contributions of the protein's N- and C-termini are accounted for similarly. 





        We recommend that the 'kernel_list' is symmetric about index 'origin' and unimodal. 
        For instance '[1,2,3,2,1]' and '2'.


        :param kernel_list: The list of how many copies of nearby residues we use when smoothing the isoelectric points
        :param origin: The index in the 'kernel_list' of the current residue
        :param inplace: Whether this modifies 'self.pI_list' or returns a new list
        :return: For 'inplace==False' the new, smoothed list of isoelectric point values. For 'inplace==True' nothing is returned.
        """


        k = kernel_list
        for a in kernel_list:
            assert type(a) == int
            assert a >= 0
        assert 0 <= origin < len(k)
        fasta_header, fasta = re.findall(string = self.fasta, pattern = r'^(>.+)\n([A-Z]*)$')[0]
        out_pI_list = []

        for i in range(len(self.pI_list)):
            local_fasta_list = []
    
            for j in range(len(k)):
                #print(i,j) #testing
                try:
                    if i+j -origin <0:
                        continue
                    local_fasta_list += [fasta[i+j -origin ] ]* k[j]
                except IndexError:
                    pass
                    
                if i+j -origin == 0:
                    N_term_count = k[j]
                    
                else:
                    N_term_count = 0
                if i+j -origin == len(self.pI_list) -1:
                    C_term_count = k[j]
                    
                else:
                    C_term_count = 0
                
            out_pI_list.append(pI_iter_alg( local_fasta_list, N_term_count = N_term_count, C_term_count = C_term_count))

        if inplace:
            self.pI_list = out_pI_list
            return 0
        else:
            return out_pI_list
            
    @staticmethod
    def run_ssearch_indices(p1: 'FGW_protein',
    p2: 'FGW_protein',
    allow_mismatch: bool = True): 
        """
        Runs the ssearch36 program from the fasta36 packages and returns the indices of the two proteins which are aligned.
        :param p1: First protein
        :param p2: Fecond protein
        :param allow_mismatch: Whether to include residues which are aligned but not the same type of amino acid
        :return: Two lists of indices, those of 'p1' and 'p2' which are aligned
        """ 

        return run_fasta36.run_ssearch_cigar_Ram(fasta1 = p1.fasta, fasta2 = p2.fasta, allow_mismatch = allow_mismatch)

    def scale_ipdm(self,
        scaler: Callable[[float],float] = lambda x :x, 
        inplace: bool = False):
        """
        :param scaler: A function with which to scale the intraprotein distance matrix. It must send 0 to 0, be strictly monotonic increasing, and concave down.
        :param inplace: Whether to modify 'self.ipdm' or output the scaled ipdm.
        :return: The scaled ipdm if 'inplace == False', and 'None' if 'inplace == True'.
        """

        m= np.vectorize(scaler)(self.ipdm)
        if inplace:
            self.ipdm = m
            self.scaled_flag = True
        else:
            return FGW_protein(fasta = self.fasta, pI_list = self.pI_list, ipdm = m, coords = None, name = self.name+'_scaled', scaled_flag = True)


    def make_GW_cell(self, 
        distribution: npt.NDArray[np.float_] = None) -> gw_cython.GW_cell:
        """
        Makes a 'gw_cython-GW_cell' object from the CAJAL softward package. This allows more efficient GW computations, but is not suitable for FGW.
        :param distribution: The mass distribution to be used.
        :return: A 'gw_cython-GW_cell' object representing 'self'
        """

        if distribution is None:
            distribution = GW_scripts.unif(len(self.pI_list))

        assert len(distribution) == self.ipdm.shape[0]
        #assert np.sum(distribution) == 1 #possibly off do to floating point rounding


        return gw_cython.GW_cell(self.ipdm,  distribution)
            

    @staticmethod       
    def run_GW_from_cells(
        cell_1:gw_cython.GW_cell  , 
        cell_2: gw_cython.GW_cell) -> float:

        """
        This is a wrapper for the CAJAL code to compute the GW distance between 'cell_1' and 'cell_2'
        :param cell_1:
        :param cell_2:
        :return: Returns the GW distance
        """

        return GW_scripts.GW_identity_init(cell_1, cell_2)
        
    @staticmethod       
    def run_GW(P1 :'FGW_protein',
        P2: 'FGW_protein') -> float:
        """
        This is a wrapper for the CAJAL code to create gw_cython.GW_cell object then compute the GW distance between them
        :param P1:
        :param P2:
        :return: Returns the GW distance
        """

        cell_1 = P1.make_GW_cell()
        cell_2 = P2.make_GW_cell()
        return run_GW_from_cells(cell_1, cell_2)


    def downsample_by_indices(self, 
        indices: list[int]) -> 'FGW_protein':
        """
        This creates a new 'FGW_protein' object consisting of the residues of 'self' in the input indices
        :param indices: The indices to keep.
        :return: A new 'FGW_protein' object.
        """
        assert set(indices).issubset(set(range(len(self.pI_list))))
        if self.coords is not None:
            coords = self.coords[indices]
        else:
            coords = None
        ii = np.ix_(indices,indices)
        ipdm = self.ipdm[ii]
        pI_list = [self.pI_list[i] for i in indices]

        fasta_header, fasta_seq = re.findall( string = self.fasta, pattern = r'^(>.+)\n([A-Z]*)$')[0]
        new_header = fasta_header + ' downsampled'
        new_seq = ''.join([fasta_seq[i] for i in indices])
        new_fasta = new_header + '\n' + new_seq
        
        
        #fasta_seq = self.fasta_seq[indices] #deprecated/unnecesary i think
        return FGW_protein(fasta = new_fasta, pI_list = pI_list, ipdm = ipdm, coords = coords, name = self.name+'_downsampled', scaled_flag = self.scaled_flag)

    def recompute_ipdm(self) -> None:
        """
        This method recalculates the ipdm based on the coordinates. The two might not be compatible because of scaling or 'downsample_n(mean_sample =True)'.
        Raises an error if 'self.coords is None'.
        :return: Does not return
        """
        if self.coords is None:
            raise Exception('self.coords is None')
        else:
            self.ipdm = squareform(pdist(self.coords))
            self.scaled_flag = False
        
    def downsample_n(self,
        n:int = np.inf, 
        pI_combination: bool = True,
        pI_alg: str = 'iter',
        left_sample:bool = False,
        mean_sample:bool = False) -> 'FGW_protein':
        """
        This method makes a new 'FGW_protein' object created by downsampling from 'self'. This is done by dividing 'self' into 'n' evenly sized segments, 
        then creates an 'FGW_protein' object whose residues are formed by those segments. Depending on the parameters this can be done with regular downsampling 
        (simply picking one residue from each segment and copying its data) or by combining the coordinate data and/or isoelectric values of the residues in a segment.

        :param n: The maximum number of residues in the output protein. If this is larger than the size of 'self', then there is no downsampling.
        :param pI_combination: Whether to combine the isoelectric points of nearby residues when downsampling. If "False" then the values in 'pI_list' 
            of the returned 'FGW_protein' are a subset of those of 'self.pI_list'. If "True" then 'pI_alg' is used to estimate the isoelectric point of
            nearby residues
        :param pI_alg: Which algorithm to use to estimate isoelectric points. 'pI_alg == "iter"' uses 'pI_iter_alg()' and 'pI_alg == "median"' uses 'pH_median'
        :param left_sample: Whether to use the left-most (lowest index) or median residue from each segment. 'left_sample == True' uses the left-most, 
            'left_sample== False' uses the median
        :param mean_sample: Whether to average the coordinates of the residues in a segment. 'mean_sample == False' uses the coordinates of the residue determined by 'left_sample',
            'mean_sample==True' uses the average of the coordinates in a segment.
        :return: A new 'FGW_protein' object created by downsampling from 'self'.
        """


        n = min(n, len(self.pI_list))
        l,s = np.linspace(0, len(self.pI_list), num=n, endpoint=False, dtype=int,retstep = True) 
        if left_sample:
            indices = np.array([int(i ) for i in l])
        else:
            indices = np.array([int(i + s//2) for i in l])


        if self.coords is not None:
            if not mean_sample:
                coords = self.coords[indices, :] #untested
    
            else: 
                split_coord_list = read_pdb.split_list(self.coords, n)
                coords = [np.mean(seg, axis = 0) for seg in split_coord_list]  #unsure about axis
                coords = np.stack(coords) 
            ipdm = None
        else:
            coords = None

        
        ii = np.ix_(indices,indices)
        ipdm = self.ipdm[ii]

        if pI_combination: 
            split_pI_list = read_pdb.split_list(self.pI_list, n)
            if pI_alg == 'median':
                pI_list = [pH_median(seg) for seg in split_pI_list]
            elif pI_alg == 'iter':
                split_res_list = read_pdb.split_list(self.get_fasta_seq(), n)
                pI_list = [pI_iter_alg(split_res_list[0], N_term_count= 1)]  + [pI_iter_alg(seg) for seg in split_res_list[1:-1]] + [pI_iter_alg(split_res_list[-1], C_term_count= 1)]
            else:
                raise Exception("Invalid parameter for pI_alg, must be 'median' or 'iter'" )
        else:
            pI_list = [self.pI_list[i] for i in indices]

        fasta_header, fasta_seq = re.findall( string = self.fasta, pattern = r'^(>.+)\n([A-Z]*)$')[0]
        new_header = fasta_header + ' downsampled to ' + str(n)
        new_seq = ''.join([fasta_seq[i] for i in indices])
        new_fasta = new_header + '\n' + new_seq
        
        
        #fasta_seq = self.fasta_seq[indices] #deprecated/unnecesary i think
        return FGW_protein(fasta = new_fasta, pI_list = pI_list, ipdm = ipdm, coords = coords, name = self.name+'_downsampled', scaled_flag = self.scaled_flag)


    def validate(self) -> bool:
        """
        Checks if a 'FGW_protein' object passes basic tests.
        :return: 'True' is it passes, raises assertion error otherwise.
        """
        if not self.coords is None:
            
            assert type(self.coords) == np.ndarray #
            assert len(self.coords.shape) == 2
            assert self.coords.shape[1] == 3 
            assert self.coords.shape[0] == len(self.pI_list) 
            if not self.scaled_flag:
                assert (self.ipdm ==squareform(pdist(self.coords))).all()

        
        assert type(self.ipdm) == np.ndarray #
        assert len(self.ipdm.shape) ==2
        assert self.ipdm.shape[0] == self.ipdm.shape[1]
        assert (self.ipdm == self.ipdm.T).all()
        assert self.name is not None

        fasta_header, fasta_seq = re.findall(string = self.fasta, pattern = r'^(>.+)\n([A-Z]*)$')[0]
        assert len(fasta_seq) ==len(self.pI_list)

        if  len(self.pI_list)>=1 and ( self.pI_list[1:-1] != [read_pdb.writeProtIepMedian(r) for r in fasta_seq[1:-1]]):
            print('pI_list is wrong, could be caused by convolution')
        
        return True

    def get_fasta_seq(self)->str:
        """
        Helper method to get the sequence of a protein
        :return: The protein sequence
        """
        fasta_header, fasta_seq = re.findall( string = self.fasta, pattern = r'^(>.+)\n([A-Z]*)$')[0]
        return fasta_seq
        
    @staticmethod
    def make_protein_from_files(pdb_file: str, 
        fasta_file:str) -> 'FGW_protein':
        """
        Creates a FGW_protein object with the coordinate data from the 'pdb_file' and the sequence of the 'fasta_file'.
        :param pdb_file: Filepath to the .pdb file
        :param fasta: Filepath to the .fasta file
        :return: A new 'FGW_protein' object
        """

        with open(fasta_file, 'r') as fasta_in:
            fasta = fasta_in.read()
        coords, pI_list = read_pdb.get_pdb_coords_pI(filepath = pdb_file, n = np.inf, median = True)
        name = re.findall(string = pdb_file, pattern = r'([^\/]+)\.pdb$')[0]
        return FGW_protein(name = name, coords = coords, pI_list = pI_list,fasta=fasta)
        
    @staticmethod
    def make_protein_from_pdb(pdb_file:str) ->'FGW_protein':
        """
        Creates a FGW_protein object with the coordinate and sequence data from the 'pdb_file'
        :param pdb_file: Filepath to the .pdb file
        :return: A new 'FGW_protein' object
        """

        coords, pI_list = read_pdb.get_pdb_coords_pI(filepath = pdb_file, n = np.inf, median = True)
        name = re.findall(string = pdb_file, pattern = r'([^\/]+)\.pdb$')[0]
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_file)
    
        # Extract sequence from structure
        sequence = ""
        for model in structure:
            for chain in model:
                for residue in chain:
                    if PDB.is_aa(residue.get_resname(), standard=True):
                        sequence += PDB.Polypeptide.protein_letters_3to1[residue.get_resname()]
                    elif 'UNK' in residue.get_resname():
                        sequence += 'X' #unkown
        assert len(sequence) == len(pI_list)
        fasta = ">" + name + '\n' + sequence    
        return FGW_protein(name = name, coords = coords, pI_list = pI_list,fasta=fasta)


    @staticmethod
    def run_FGW(p1: 'FGW_protein', p2:'FGW_protein', alpha:float) -> float:
        """
        This calculates the fused Gromov-Wasserstein distance between two proteins. The computation is done with the Python 'ot' library. 
        :param p1: The first protein
        :param p2: The second protein
        :param alpha: The trade-off parameter in [0,1] between fused term and geometric term. A higher value of 'alpha' means more geometric weight, 'alpha' = 1 is equivalent to regular GW.
        :return: The FGW distance
        """

        D1 = p1.ipdm
        D2 = p2.ipdm
        pI1 = p1.pI_list
        pI2 = p2.pI_list
        n1 = len(D1)
        n2 = len(D2)
        try:
            assert n1 == len(pI1)
            assert n2 == len(pI2)
        except:
            print(D1.shape, D2.shape, len(pI1), len(pI2))
            assert False
        
        a = np.array([np.array([x]) for x in pI1])
        b = np.array(pI2)
        aa = np.broadcast_to(a,(n1,n2))
        bb = np.broadcast_to(b,(n1,n2))
        M = abs(aa-bb)
        G0 = GW_scripts.id_initial_coupling_unif(n1,n2)
        
        d = ot.fused_gromov_wasserstein2(M=M, C1=D1, C2=D2, alpha = alpha, p= GW_scripts.unif(n1),q=GW_scripts.unif(n2), G0 = G0, loss_fun='square_loss')
    
    
        return  0.5 * math.sqrt(d)
        
    @staticmethod
    def run_FGW_seq_aln(p1:'FGW_protein', p2:'FGW_protein', alpha:float, n: int = np.inf,allow_mismatch:bool = True) -> float:
        """
        This calculates the fused Gromov-Wasserstein distance between two proteins when applied just to aligned residues. 
        It first applies sequence alignment, downsamples up to 'n' of the aligned residues, then applies FGW. 
        :param p1: The first protein
        :param p2: The second protein
        :param n: The maximum number of residues to use (to reduce runtime).
        :param alpha: The trade-off parameter in [0,1] between fused term and geometric term. A higher value of 'alpha' means more geometric weight, 'alpha' = 1 is equivalent to regular GW.
        :return: The FGW distance
        """
        inds1, inds2 = FGW_protein.run_ssearch_indices(p1 =p1, p2 = p2, allow_mismatch = allow_mismatch)

        if n < len(inds1):
            l,s = np.linspace(0, len(inds1), num=n, endpoint=False, dtype=int,retstep = True) #new method I'm trying
            subindices = np.array([int(i + s//2) for i in l])

            inds1 = [inds1[i] for i in subindices]
            inds2 = [inds2[i] for i in subindices]
            
        p3 = p1.downsample_by_indices(inds1)
        p4 = p2.downsample_by_indices(inds2)

        return FGW_protein.run_FGW(p3,p4, alpha = alpha)
        




    

    