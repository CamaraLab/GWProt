import os
#os.environ['OPENBLAS_NUM_THREADS'] = '1' #should be first, before ot

from threadpoolctl import ThreadpoolController
controller = ThreadpoolController()
import re
import math
import numpy as np
import ot
import numpy.typing as npt
from scipy.spatial.distance import *

import Bio.PDB
from Bio import PDB, SeqIO
from typing import Callable, Union



from .GW_scripts import *
from .read_pdb import *
from .run_fasta36 import *

from cajal import run_gw, qgw, gw_cython
"""
copied 4/1/2024

"""



class GW_protein:
    """
    This class contains everything needed to run GW and FGW on proteins, as well as versions with distortion scaling and sequence alignment
    
    :param name: A string for ease of use
    :param coords: The coordinates of the CA atoms of the protein, ordered sequentially
    :param seq: A string giving the sequence of the protein
    :param ipdm: The intra-protein distance matrix of a protein. The (i,j)th entry is the (possibly scaled) distance between residues i and j. This is mutable can can change if distortion scaling is used.
    :param scaled_flag: Records whether the ipdm is the exact distance between residues or if it has been scaled.
    :param distribution: np.array of the weighting of the residues, must sum to 1. Default is a uniform distribution.

    """




    def __init__(self, name : str, seq: str,  coords = None, ipdm = None, scaled_flag :bool = False , distribution = None):

        assert not (coords is None and ipdm is None)
        if not coords is None:
            coords = np.array(coords)


            
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
        Compares the ``seq``, the ``ipdm``, ``distribution``, and the ``coords`` if both are defined.
        This does *NOT* compare the ``name`` or ``scaled_flag``.
        """
        
        if self.coords is not None and other.coords is not None and ((self.coords.shape != other.coords.shape) or (self.coords != other.coords).any()):
            return False  
        return self.seq == other.seq  and (self.ipdm == other.ipdm).all() and (self.distribution == other.distribution).all()
      

    def __len__(self):
        """:return: the number of amino acids in the protein"""
        return self.ipdm.shape[0]
        
    def __str__(self):
        return self.name
 
    @controller.wrap(limits=1, user_api=controller.info()[-1]['user_api'])
    @staticmethod
    def run_ssearch_indices(prot1: 'GW_protein',
    prot2: 'GW_protein',
    allow_mismatch: bool = True) -> tuple[list[int],list[int]]: 

        """

        Runs a local sequence alignment returns the indices of the two proteins which are aligned. ssearch36 must be in the PATH to use this method.
        
        :param prot1: First protein
        :param prot2: Fecond protein
        :param allow_mismatch: Whether to include residues which are aligned but not the same type of amino acid
        :return: Two lists of indices, those of ``prot1`` and ``prot2`` which are aligned

        """ 
        fasta1 = ">" + prot1.name + '\n' + prot1.seq
        fasta2 = ">" + prot2.name + '\n' + prot2.seq
        return run_ssearch_cigar_Ram(fasta1 = fasta1, fasta2 = fasta2, allow_mismatch = allow_mismatch)

    def scale_ipdm(self,
        scaler: Callable[[float],float] = math.sqrt, 
        inplace: bool = False):

        """
        This method scales all entries of the intra-protein distance matrix.

        :param scaler: A function with which to scale the intraprotein distance matrix. It must send 0 to 0, be strictly monotonic increasing, and concave down. 
            Default is the square root function.
        :param inplace: Whether to modify ``self.ipdm`` or output a new ``GW_protein`` object.
        :return: The scaled ipdm if ``inplace == False``, and ``None`` if ``inplace == True``.

        """

        m= np.vectorize(scaler)(self.ipdm)
        if inplace:
            self.ipdm = m
            self.scaled_flag = True
        else:
            return GW_protein(seq = self.seq, ipdm = m, coords = None, name = self.name+'_scaled', scaled_flag = True, distribution = self.distribution)


    def make_cajal_cell(self) -> gw_cython.GW_cell:

        """
        This method makes a ``cajal.gw_cython.GW_cell`` object from the CAJAL library. 

        :return: A ``cajal.gw_cython.GW_cell`` object representing ``self``.

        """

        return gw_cython.GW_cell(self.ipdm,  self.distribution)
            
    @controller.wrap(limits=1, user_api=controller.info()[-1]['user_api'])
    @staticmethod       
    def run_GW_from_cajal(
        cajal_cell1:gw_cython.GW_cell  , 
        cajal_cell2: gw_cython.GW_cell,
        transport_plan:bool = False) -> Union[float, tuple[float, np.array]]:

        """
        This is a wrapper for the CAJAL code to compute the GW distance between ``cajal_cell1`` and ``cajal_cell2``, 
        outputs the computed transport plan if ``tranport_plan``.

        :param cajal_cell1:
        :param cajal_cell2:
        :param transport_plan: Whether to return the computed transport plan
        :return: Returns the GW distance and optimal transport plan if ``transport_plan``

        """

        return GW_identity_init(cajal_cell1, cajal_cell2, transport_plan= transport_plan)
        
    @controller.wrap(limits=1, user_api=controller.info()[-1]['user_api'])
    @staticmethod       
    def run_GW(prot1 :'GW_protein',
        prot2: 'GW_protein',
        transport_plan:bool = False) -> Union[float, tuple[float, np.array]]:

        """
        Computes the GW distance and transport plan if ``transport_plan``.

        :param prot1:
        :param prot2:
        :param transport_plan: Whether to return the computed transport plan
        :return: Returns the GW distance and optimal transport plan if ``transport_plan``
        """

        P1, P2 = prot1, prot2

        cell_1 = P1.make_cajal_cell()
        cell_2 = P2.make_cajal_cell()
        return GW_protein.run_GW_from_cajal(cell_1, cell_2, transport_plan = transport_plan)


    def downsample_by_indices(self, 
        indices: list[int]) -> 'GW_protein':

        """
        This creates a new ``GW_protein`` object consisting of the residues of ``self`` in the input indices.

        :param indices: The indices to keep.
        :return: A new ``GW_protein`` object
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
        return GW_protein(seq = new_seq,  ipdm = ipdm, coords = coords, name = self.name+'_downsampled', scaled_flag = self.scaled_flag, distribution = new_distribution)

    def reset_ipdm(self) -> None:

        """
        This method recalculates the ipdm inplace based on the coordinates. 
        Raises an error if ``self.coords is None``.
        """

        if self.coords is None:
            raise Exception('self.coords is None')
        else:
            self.ipdm = squareform(pdist(self.coords))
            self.scaled_flag = False
        #The two might not be compatible because of scaling or ``downsample_n(mean_sample =True)``.
        
  

    def validate(self) -> bool:

        """
        Checks if a ``GW_protein`` object passes basic consistency tests.

        :return: ``True`` is it passes, raises assertion error otherwise.
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
    def make_protein_from_pdb(pdb_file:str, chain_id: str = None) ->'GW_protein':

        """
        Creates a ``GW_protein`` object with the coordinate and sequence data from the ``pdb_file``. This gives a uniform distribution.
        
        :param pdb_file: Filepath to the pdb file
        :param chain_id: Which chain(s) to use, None uses all chains
        :return: A new ``GW_protein`` object

        """

        coords, _ , seq = get_pdb_coords_pI(filepath = pdb_file, n = np.inf, median = True, chain_id = chain_id)
        name = re.findall(string = pdb_file, pattern = r'([^\/]+)\.pdb$')[0]
        if chain_id is not None:
            name += '_'+chain_id

        
        return GW_protein(name = name, coords = coords, seq=seq)


    @controller.wrap(limits=1, user_api=controller.info()[-1]['user_api'])
    @staticmethod
    def run_FGW_data_lists(prot1: 'GW_protein', prot2:'GW_protein', data1 :list[float] , data2 : list[float] , alpha:float = 1, transport_plan = False) -> Union[float, tuple[float, np.array]]:
        
        """
        This calculates the fused Gromov-Wasserstein distance between two proteins. 
        It takes in a list of ``float`` s for each proteins representing the value in the feature space for each residue. 
        The ijth entry in the associated distance matrix is ``abs(data1[i] - data2[j])``.
        
        :param prot1: The first protein
        :param prot2: The second protein
        :param data1: The data used in the first protein
        :param data2: The data used in the second protein
        :param alpha: The trade-off parameter in [0,1] between fused term and geometric term. A higher value of ``alpha`` means more geometric weight, 
            ``alpha=1`` is equivalent to regular GW.
        :param transport_plan: Whether to return the computed transport plan
        :return: Returns the FGW distance and transport plan if ``transport_plan``

        """

        p1,p2 = prot1, prot2
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
        G0 = id_initial_coupling(p1.distribution,p2.distribution)

        T , log= ot.fused_gromov_wasserstein(M=M, C1=D1, C2=D2, alpha = alpha, p= p1.distribution ,q=p2.distribution, G0 = G0, loss_fun='square_loss', log = True)
        d = 0.5 * math.sqrt(max(0,log['fgw_dist']))

        if transport_plan:
            return d, T
        else:
            return d
            

    def downsample_n(self,
        n:int = np.inf, 
        left_sample:bool = False,
        mean_sample:bool = False) -> 'GW_protein':

        """
        This method makes a new ``GW_protein`` object created by downsampling from ``self``. This is done by dividing ``self`` into ``n`` evenly sized segments, 
        then creates an ``GW_protein`` object whose residues are formed by those segments.

        :param n: The maximum number of residues in the output protein. If this is larger than ``len(self)``, 
            then there is no downsampling.
        :param left_sample: Whether to use the left-most (lowest index) or median residue from each segment. ``left_sample == True`` uses the left-most, ``left_sample== False`` uses the median.
        :param mean_sample: Whether to average the coordinates of the residues in a segment. ``mean_sample == False`` uses the coordinates of the residue determined by ``left_sample``, ``mean_sample==True`` uses the average of the coordinates in a segment.
        :return: A new ``GW_protein`` object created by downsampling from ``self``.

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
                split_coord_list = split_list(self.coords, n)
                coords = [np.mean(seg, axis = 0) for seg in split_coord_list]  #unsure about axis
                coords = np.stack(coords) 
            ipdm = None
        else:
            coords = None

        
        ii = np.ix_(indices,indices)
        ipdm = self.ipdm[ii]

        


        new_seq = ''.join([self.seq[i] for i in indices])

        new_distribution = np.array([ np.sum(seg) for seg in  split_list(list(self.distribution), n) ])
        
        
        return GW_protein(seq = new_seq, ipdm = ipdm, coords = coords, name = self.name+'_downsampled', scaled_flag = self.scaled_flag, distribution = new_distribution)



    
    def get_eccentricity(self, p: float =2) ->np.array:

        """
        This calculates the eccentricity of residues in a protein with exponent ``p``.
        
        :param p: The exponent, ``0< p <= np.inf``
        :return: The eccentricities of each residue, as a ``np.array``

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
    def GW_stress(prot1: 'GW_protein', prot2: 'GW_protein', T: np.array) -> tuple[np.array,np.array]:

        """
        This calculates the stress, i.e. the contribution of each residue to the sum in the GW cost, using the transport plan ``T``.
        This is output as two ``np.array`` s, one for ``prot1`` , the second for ``prot2``.

        :param prot1: The first ``GW_protein``
        :param prot2: The second ``GW_protein``
        :param T: The transport plan to be used
        :return: ``stress1, stress2``; the stresses for the two proteins

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
    def FGW_stress(prot1: 'GW_protein', prot2: 'GW_protein', T: np.array, diff_mat : np.array, alpha:float)-> tuple[np.array,np.array]:
        
        """
        This calculates the stress, i.e. the contribution of each residue to the sum in the FGW cost, using the transport plan ``T``.
        This is output as two ``np.array`` s, one for ``prot1`` , the second for ``prot2``.

        :param prot1: The first ``GW_protein``
        :param prot2: The second ``GW_protein``
        :param diff_mat: The difference matrix in the feature space
        :param T: The transport plan to be used
        :param alpha: The trade-off constant between the fused cost and the geometric cost
        :return: ``stress1, stress2``;  the stresses for the two proteins

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
        
    @controller.wrap(limits=1, user_api=controller.info()[-1]['user_api'])
    @staticmethod
    def run_FGW_dict(prot1: 'GW_protein', prot2:'GW_protein', d: dict[str , dict[str ,float]] ,alpha:float = 1, transport_plan: bool = False) -> Union[float, tuple[float, np.array]]:
        
        """
        This calculates the fused Gromov-Wasserstein distance between two proteins.  
        
        :param prot1: The first protein
        :param prot2: The second protein
        :param d: The dictionary used for the fused distances based on the protein sequences. Of the form ``d['A']['B'] == float``
        :param alpha: The trade-off parameter in [0,1] between fused term and geometric term. A higher value of ``alpha`` means more geometric weight, ``alpha = 1`` is equivalent to regular GW.
        :param transport_plan: Whether to return the transport plan
        :return: Returns the FGW distance and transport plan if ``transport_plan``

        """

        assert 0 <= alpha <=1
        D1 = prot1.ipdm
        D2 = prot2.ipdm
        n1 = len(prot1)    
        n2 = len(prot2)
        # M has shape (n1,n2)
        M = np.zeros((n1,n2))
        for i in range(n1):
            for j in range(n2):
                if prot1.seq[i] in d.keys() and prot2.seq[j] in d[prot1.seq[i]].keys():
                    M[i,j] = d[prot1.seq[i]][prot2.seq[j]]
                else:
                    M[i,j] = np.nan

        #must be tested:
        col_means = np.nanmean(M, axis = 0)
        inds = np.where(np.isnan(M))
        M[inds] = np.take(col_means, inds[1])
        row_means = np.nanmean(M, axis = 1)
        inds = np.where(np.isnan(M))
        M[inds] = np.take(row_means, inds[0])
        full_mean = np.nanmean(M)
        inds = np.where(np.isnan(M))
        M[inds] = full_mean


        G0 = id_initial_coupling(prot1.distribution,prot2.distribution)

        T , log= ot.fused_gromov_wasserstein(M=M, C1=D1, C2=D2, alpha = alpha, p= prot1.distribution ,q=prot2.distribution, G0 = G0, loss_fun='square_loss', log = True)
        d = 0.5 * math.sqrt(max(0,log['fgw_dist']))

        if transport_plan:
            return d, T
        else:
            return d


    @controller.wrap(limits=1, user_api=controller.info()[-1]['user_api'])
    @staticmethod
    def run_FGW_diff_mat(prot1: 'GW_protein', prot2:'GW_protein', diff_mat: np.array ,alpha:float = 1, transport_plan: bool = False) -> Union[float, tuple[float, np.array]]:
        
        """
        This calculates the fused Gromov-Wasserstein distance between two proteins. 

        :param prot1: The first protein
        :param prot2: The second protein
        :param diff_mat: A user-inputted matrix of the differences in the feature space between the residues of the two proteins. Of shape ``(len(prot1),len(prot2))``. \
            ``diff_mat[i,j]`` is the difference in features between the ith residue of ``prot1`` and the jth residue of ``prot2``.
        :param alpha: The trade-off parameter in [0,1] between fused term and geometric term. 
            A higher value of ``alpha`` means more geometric weight, ``alpha = 1`` is equivalent to regular GW.
        :param transport_plan: Whether to return the transport plan
        :return: Returns the FGW distance and optimal transport plan if ``transport_plan``

        """
        p1,p2 = prot1, prot2
        assert 0 <= alpha <=1
        D1 = p1.ipdm
        D2 = p2.ipdm
        n1 = len(p1)
        n2 = len(p2)
        assert diff_mat.shape == (n1,n2)

        
        G0 = id_initial_coupling(p1.distribution,p2.distribution)

        T , log= ot.fused_gromov_wasserstein(M=M, C1=D1, C2=D2, alpha = alpha, p= p1.distribution ,q=p2.distribution, G0 = G0, loss_fun='square_loss', log = True)
        d = 0.5 * math.sqrt(max(0,log['fgw_dist']))

        if transport_plan:
            return d, T
        else:
            return d
    

        


    @controller.wrap(limits=1, user_api=controller.info()[-1]['user_api'])
    @staticmethod
    def run_GW_seq_aln(prot1:'GW_protein', prot2:'GW_protein',  allow_mismatch:bool = True) -> float:
        """
        This calculates the Gromov-Wasserstein distance between two proteins when applied just to aligned residues. 
        It first applies sequence alignment, downsamples to the aligned residues, then applies GW. ssearch36 must be in the PATH to use this method.

        :param prot1: The first protein
        :param prot2: The second protein
        :return: The GW distance

        """

        p1, p2 = prot1 , prot2

        inds1, inds2 = GW_protein.run_ssearch_indices(p1 =p1, p2 = p2, allow_mismatch = allow_mismatch)


            
        p3 = p1.downsample_by_indices(inds1)
        p4 = p2.downsample_by_indices(inds2)

        return GW_protein.run_FGW(p3,p4, transport_plan = False)

    