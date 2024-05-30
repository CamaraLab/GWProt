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
        #note - the fasta is the sequence, not the file
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
      

    def __len__(self):
        return len(self.pI_list)
 
            
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
                        sequence += '*' #unkown
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
        





    