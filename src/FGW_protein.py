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
import sparse

import Bio.PDB
from Bio import PDB, SeqIO
from typing import Callable



import GW_scripts
import read_pdb
import run_fasta36
import FGW_matrices
import pymol_protein_viewer

from cajal import run_gw, qgw, gw_cython
"""
copied 4/1/2024

"""



class FGW_protein:
    """
    This class contains everything needed to run GW and fused GW on proteins, as well as versions with distortion scaling and sequence alignment

    :param name: Simply for ease of use
    :param coords: The coordinates of the CA atoms of the protein, ordered sequentially
    :param seq: A string giving the sequence of the protein
    :param ipdm: The intra-protein distance matrix of a protein. 
        The (i,j)th entry is the (possibly scaled) distance between residues i and j. This is mutable can can change if distortion scaling is used.
    :param scaled_flag: Records whether the ipdm is the exact distance between residues or if it has been scaled.
    :param distribution: Numpy array of the weighting of the residues, must sum to 1. Default is a uniform distribution.

    """




    def __init__(self, name : str, seq: str,  coords = None, ipdm = None, scaled_flag :bool = False , distribution = None):

        assert not (coords is None and ipdm is None)
        if not coords is None:
            coords = np.array(coords)

            # print('type(coords)' ,type(coords)) #debugging
            # print('coords.shape' ,coords.shape)
            
            assert len(coords.shape) == 2
            assert coords.shape[1] == 3 
        
        if ipdm is not None:
            ipdm = np.array(ipdm)
            assert len(ipdm.shape) ==2
            assert ipdm.shape[0] == ipdm.shape[1]
            assert (ipdm == ipdm.T).all()

        
        self.name = name
        self.seq = seq
        
        self.coords = coords


        
        self.scaled_flag = scaled_flag #whether the ipdm has been scaled

        if ipdm is None:
            self.ipdm = squareform(pdist(self.coords))
        else:
            
            self.ipdm = ipdm

        if distribution is None:
            self.distribution = np.ones(self.ipdm.shape[0])/ self.ipdm.shape[0]
        else:
            assert distribution.shape[0] == self.ipdm.shape[0]
            self.distribution = distribution
        assert math.isclose(np.sum(self.distribution),1)
            
    def __eq__(self, other):
        """
        Compares the underlying sequences, the ipdms, distributions, and the coords if both are defined.
        This does NOT compare the names or scaled_flags.
        """
        
        if self.coords is not None and other.coords is not None and ((self.coords.shape != other.coords.shape) or (self.coords != other.coords).any()):
            return False  
        return self.seq == other.seq  and (self.ipdm == other.ipdm).all() and (self.distribution == other.distribution).all()
      

    def __len__(self):
        return self.ipdm.shape[0]
        
    def __str__(self):
        return self.name
 
            
    @staticmethod
    def run_ssearch_indices(p1: 'FGW_protein',
    p2: 'FGW_protein',
    allow_mismatch: bool = True): 
        """
        Runs the ssearch36 program from the seq36 packages and returns the indices of the two proteins which are aligned.
        :param p1: First protein
        :param p2: Fecond protein
        :param allow_mismatch: Whether to include residues which are aligned but not the same type of amino acid
        :return: Two lists of indices, those of 'p1' and 'p2' which are aligned
        """ 
        fasta1 = ">" + p1.name + '\n' + p1.seq
        fasta2 = ">" + p2.name + '\n' + p2.seq
        return run_fasta36.run_ssearch_cigar_Ram(fasta1 = fasta1, fasta2 = fasta2, allow_mismatch = allow_mismatch)

    def scale_ipdm(self,
        scaler: Callable[[float],float] = math.sqrt, 
        inplace: bool = False):
        """
        :param scaler: A function with which to scale the intraprotein distance matrix. It must send 0 to 0, be strictly monotonic increasing, and concave down. 
        Default is the square root function.
        :param inplace: Whether to modify 'self.ipdm' or output the scaled ipdm.
        :return: The scaled ipdm if 'inplace == False', and 'None' if 'inplace == True'.
        """

        m= np.vectorize(scaler)(self.ipdm)
        if inplace:
            self.ipdm = m
            self.scaled_flag = True
        else:
            return FGW_protein(seq = self.seq, ipdm = m, coords = None, name = self.name+'_scaled', scaled_flag = True, distribution = self.distribution)


    def make_GW_cell(self, 
        distribution: npt.NDArray[np.float_] = None) -> gw_cython.GW_cell:
        """
        Makes a 'gw_cython-GW_cell' object from the CAJAL softward package. This allows more efficient GW computations, but is not suitable for FGW.
        :param distribution: The mass distribution to be used.
        :return: A 'gw_cython-GW_cell' object representing 'self'
        """

        return gw_cython.GW_cell(self.ipdm,  self.distribution)
            

    @staticmethod       
    def run_GW_from_cells(
        cell_1:gw_cython.GW_cell  , 
        cell_2: gw_cython.GW_cell,
        transport_plan:bool = False) -> float:

        """
        This is a wrapper for the CAJAL code to compute the GW distance between 'cell_1' and 'cell_2', 
        outputs the computed transport plan if 'tranport_plan'.
        :param cell_1:
        :param cell_2:
        :param transport_plan: Whether to return the computed transport plan
        :return: Returns the GW distance and transport plan if 'transport_plan'
        """

        return GW_scripts.GW_identity_init(cell_1, cell_2, transport_plan= transport_plan)
        
    @staticmethod       
    def run_GW(P1 :'FGW_protein',
        P2: 'FGW_protein',
        transport_plan:bool = False) -> float:
        """
        This is a wrapper for the CAJAL code to create gw_cython.GW_cell objects then compute the GW distance between them.
        Returns the GW distance and transport plan if 'transport_plan'
        :param P1:
        :param P2:
        :param transport_plan: Whether to return the computed transport plan
        :return: Returns the GW distance and transport plan if 'transport_plan'
        """

        cell_1 = P1.make_GW_cell()
        cell_2 = P2.make_GW_cell()
        return FGW_protein.run_GW_from_cells(cell_1, cell_2, transport_plan = transport_plan)


    def downsample_by_indices(self, 
        indices: list[int]) -> 'FGW_protein':
        """
        This creates a new 'FGW_protein' object consisting of the residues of 'self' in the input indices
        :param indices: The indices to keep.
        :return: A new 'FGW_protein' object
        """
        assert set(indices).issubset(set(range(self.ipdm.shape[0])))
        if self.coords is not None:
            coords = self.coords[indices]
        else:
            coords = None
        ii = np.ix_(indices,indices)
        ipdm = self.ipdm[ii]


        new_seq = ''.join([self.seq[i] for i in indices])
        new_distribution = self.distribution[indices]/np.sum(self.distribution[indices])
        return FGW_protein(seq = new_seq,  ipdm = ipdm, coords = coords, name = self.name+'_downsampled', scaled_flag = self.scaled_flag, distribution = new_distribution)

    def reset_ipdm(self) -> None:
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
            if not self.scaled_flag:
                assert (self.ipdm ==squareform(pdist(self.coords))).all()

        
        assert type(self.ipdm) == np.ndarray #
        assert len(self.ipdm.shape) ==2
        assert self.ipdm.shape[0] == self.ipdm.shape[1]
        assert (self.ipdm == self.ipdm.T).all()
        assert self.name is not None
        assert self.distribution[0] == self.ipdm.shape[1]
        assert math.isclose(np.sum(self.distribution),1) 

  
        return True


 
    @staticmethod
    def make_protein_from_pdb(pdb_file:str, chain_id:str = None) ->'FGW_protein':
        """
        Creates a FGW_protein object with the coordinate and sequence data from the 'pdb_file'
        :param pdb_file: Filepath to the .pdb file
        :param chain_id: Which chain(s) to use, None uses all chains
        :return: A new 'FGW_protein' object

        Gives uniform distribution
        """

        coords, _ , seq = read_pdb.get_pdb_coords_pI(filepath = pdb_file, n = np.inf, median = True, chain_id = chain_id)
        name = re.findall(string = pdb_file, pattern = r'([^\/]+)\.pdb$')[0]
        if chain_id is not None:
            name += '_'+chain_id

        
        return FGW_protein(name = name, coords = coords, seq=seq)


    @staticmethod
    def run_FGW_data_lists(p1: 'FGW_protein', p2:'FGW_protein', data1 :list[float] , data2 : list[float] , alpha:float = 1, transport_plan = False) -> float:
        """
        This calculates the fused Gromov-Wasserstein distance between two proteins. The computation is done with the Python 'ot' library. 
        :param p1: The first protein
        :param p2: The second protein
        :param data1: The data used in the first protein
        :param data2: The data used in the second protein
        :param alpha: The trade-off parameter in [0,1] between fused term and geometric term. A higher value of 'alpha' means more geometric weight, 'alpha' = 1 is equivalent to regular GW.
        :param transport_plan: Whether to return the computed transport plan
        :return: Returns the FGW distance and transport plan if 'transport_plan'
        """
        #not yet tested
        D1 = p1.ipdm
        D2 = p2.ipdm

        n1 = len(D1)
        n2 = len(D2)
        try:
            assert n1 == len(data1)
            assert n2 == len(data2)
        except:
            print(D1.shape, D2.shape, len(data1), len(data2))
            assert False
        
        a = np.array([np.array([x]) for x in data1])
        b = np.array(data2)
        aa = np.broadcast_to(a,(n1,n2))
        bb = np.broadcast_to(b,(n1,n2))
        M = abs(aa-bb)
        G0 = GW_scripts.id_initial_coupling(p1.distribution,p2.distribution)

        T , log= ot.fused_gromov_wasserstein(M=M, C1=D1, C2=D2, alpha = alpha, p= p1.distribution ,q=p2.distribution, G0 = G0, loss_fun='square_loss', log = True)
        d = 0.5 * math.sqrt(log['fgw_dist'])

        if transport_plan:
            return d, T
        else:
            return d
            

    def downsample_n(self,
        n:int = np.inf, 
        left_sample:bool = False,
        mean_sample:bool = False) -> 'FGW_protein':
        """
        This method makes a new 'FGW_protein' object created by downsampling from 'self'. This is done by dividing 'self' into 'n' evenly sized segments, 
        then creates an 'FGW_protein' object whose residues are formed by those segments. Depending on the parameters this can be done with regular downsampling 
        (simply picking one residue from each segment and copying its data) or by combining the coordinate data and/or isoelectric values of the residues in a segment.

        :param n: The maximum number of residues in the output protein. If this is larger than the size of 'self', then there is no downsampling.
        :param left_sample: Whether to use the left-most (lowest index) or median residue from each segment. 'left_sample == True' uses the left-most, 
            'left_sample== False' uses the median
        :param mean_sample: Whether to average the coordinates of the residues in a segment. 'mean_sample == False' uses the coordinates of the residue determined by 'left_sample',
            'mean_sample==True' uses the average of the coordinates in a segment.
        :return: A new 'FGW_protein' object created by downsampling from 'self'.
        """


        n = min(n, len(self))
        l,s = np.linspace(0, len(self), num=n, endpoint=False, dtype=int,retstep = True) 
        if left_sample:
            indices = np.array([int(i ) for i in l])
        else:
            indices = np.array([int(i + s//2) for i in l])


        if self.coords is not None:
            if not mean_sample:
                coords = self.coords[indices, :] 
    
            else: 
                split_coord_list = read_pdb.split_list(self.coords, n)
                coords = [np.mean(seg, axis = 0) for seg in split_coord_list]  #unsure about axis
                coords = np.stack(coords) 
            ipdm = None
        else:
            coords = None

        
        ii = np.ix_(indices,indices)
        ipdm = self.ipdm[ii]

        


        new_seq = ''.join([self.seq[i] for i in indices])

        new_distribution = np.array([ np.sum(seg) for seg in  read_pdb.split_list(list(self.distribution), n) ])
        
        
        return FGW_protein(seq = new_seq, ipdm = ipdm, coords = coords, name = self.name+'_downsampled', scaled_flag = self.scaled_flag, distribution = new_distribution)



    
    def get_eccentricity(self, p: float =2):
        """
        This calculates the eccentricity of a protein with exponent p as defined in https://www.math.ucdavis.edu/~saito/data/acha.read.w12/memoli-gromov-dist.pdf, Definition 5.3.
        :param p: The exponent, 0< p <= np.inf
        :return: The eccentricities of each residue, as a np.array
        """

        ipdm = self.ipdm
        distr = self.distribution

        assert p >0
        n = ipdm.shape[0]    

            
        if p == np.inf:
            ipdm_pp =  ipdm* (distr != 0)
            eccentricity = np.max(ipdm_pp, axis = 0)
        else:
            ipdm_pp = ipdm**p
            ipdm_pp_w = ipdm_pp * distr 
            pre_stress = np.sum(ipdm_pp_w, axis = 0) 
            eccentricity = pre_stress**(1/p)

        return eccentricity


    @staticmethod
    def GW_stress(prot1: 'FGW_protein', prot2: 'FGW_protein', T: np.array):
        """
        This calculates the stress, i.e. the contribution of each residue to the sum in the GW cost, using the transport plan 'T'.
        This is output as two np.arrays, one for 'prot1', the second for 'prot2'.
        WARNING - np.sum(stress1) != c, where c is the GW cost; rather math.sqrt(np.sum(stress1))/2 == c; and similarly for stress2.

        :param prot1: The first FGW_protein
        :param prot2: The second FGW_protein
        :param T: The transport plan to be used
        :return: stress1, stress2; the stresses for the two proteins
        """

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


    @staticmethod
    def FGW_stress(prot1: 'FGW_protein', prot2: 'FGW_protein', T: np.array, diff_mat : np.array, alpha:float):
        """
        This calculates the stress, i.e. the contribution of each residue to the sum in the FGW cost, using the transport plan 'T'.
        This is output as two np.arrays, one for 'prot1', the second for 'prot2'.
        WARNING - np.sum(stress1) != c, where c is the GW cost; rather math.sqrt(np.sum(stress1))/2 == c; and similarly for stress2.

        :param prot1: The first FGW_protein
        :param prot2: The second FGW_protein
        :param diff_mat: The difference matrix in the feature space
        :param T: The transport plan to be used
        :param alpha: The trade-off constant between the fused cost and the geometric cost
        :return: stress1, stress2; the stresses for the two proteins
        """        
        n1= len(prot1)
        n2 = len(prot2)
        assert T.shape == (n1,n2)
        assert diff_mat.shape == (n1,n2)
        assert 0 <= alpha <= 1


        A = prot1.ipdm
        a = prot1.distribution
        B = prot2.ipdm
        b = prot2.distribution
        M = diff_mat



        geo_stress1 = alpha * (np.einsum('ik,il->i',T,(np.einsum('ij,ij->ij',A,A) @T) )   + T @ np.einsum('kl,kl->kl',B,B) @b  -(2 * np.einsum('ab,ab->a', A @T @B, T)))
        geo_stress2 = alpha *(np.einsum('kj,lj->j' ,T @ np.einsum('kl,kl->kl',B,B), T) +a.T @ np.einsum('ij,ij->ij',A,A) @T -(2 * np.einsum('ab,ab->b', A @T @B, T)))

        


        fused_stress1 = (1-alpha) * np.einsum('ik,ik->i', M,T)
        fused_stress2 = (1-alpha) * np.einsum('ik,ik->k', M,T)

        stress1 = geo_stress1 + fused_stress1
        stress2 = geo_stress2 + fused_stress2

        return stress1, stress2
        
    @staticmethod
    def run_FGW_dict(p1: 'FGW_protein', p2:'FGW_protein', d: dict[str , dict[str ,float]] ,alpha:float = 1, transport_plan: bool = False) -> float:
        """
        This calculates the fused Gromov-Wasserstein distance between two proteins. The computation is done with the Python 'ot' library. 
        :param p1: The first protein
        :param p2: The second protein
        :param d: The dictionary used for the fused distances based on the protein sequences. Of the form 'd['A']['B'] == float'
        :param alpha: The trade-off parameter in [0,1] between fused term and geometric term. 
        A higher value of 'alpha' means more geometric weight, 'alpha' = 1 is equivalent to regular GW.
        :param transport_plan: Whether to return the transport plan
        :return: Returns the FGW distance and transport plan if 'transport_plan'
        """

        assert 0 <= alpha <=1
        D1 = p1.ipdm
        D2 = p2.ipdm
        n1 = len(p1)
        n2 = len(p2)
        # M has shape (n1,n2)
        M = np.zeros((n1,n2))
        for i in range(n1):
            for j in range(n2):
                M[i,j] = d[p1.seq[i]][p2.seq[j]]
        
        G0 = GW_scripts.id_initial_coupling(p1.distribution,p2.distribution)

        T , log= ot.fused_gromov_wasserstein(M=M, C1=D1, C2=D2, alpha = alpha, p= p1.distribution ,q=p2.distribution, G0 = G0, loss_fun='square_loss', log = True)
        d = 0.5 * math.sqrt(log['fgw_dist'])

        if transport_plan:
            return d, T
        else:
            return d


    @staticmethod
    def run_FGW_diff_mat(p1: 'FGW_protein', p2:'FGW_protein', diff_mat: np.array ,alpha:float = 1, transport_plan: bool = False) -> float:
        """
        This calculates the fused Gromov-Wasserstein distance between two proteins. The computation is done with the Python 'ot' library. 
        :param p1: The first protein
        :param p2: The second protein
        :param d: A user-inputted matrix of the differences in the feature space between the residues of the two proteins. Of shape ('len(p1)','len(p2)').
        :param alpha: The trade-off parameter in [0,1] between fused term and geometric term. 
        A higher value of 'alpha' means more geometric weight, 'alpha' = 1 is equivalent to regular GW.
        :param transport_plan: Whether to return the transport plan
        :return: Returns the FGW distance and transport plan if 'transport_plan'
        """

        assert 0 <= alpha <=1
        D1 = p1.ipdm
        D2 = p2.ipdm
        n1 = len(p1)
        n2 = len(p2)
        assert diff_mat.shape == (n1,n2)

        
        G0 = GW_scripts.id_initial_coupling(p1.distribution,p2.distribution)

        T , log= ot.fused_gromov_wasserstein(M=M, C1=D1, C2=D2, alpha = alpha, p= p1.distribution ,q=p2.distribution, G0 = G0, loss_fun='square_loss', log = True)
        d = 0.5 * math.sqrt(log['fgw_dist'])

        if transport_plan:
            return d, T
        else:
            return d
    

        


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




    