
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1' 

from threadpoolctl import ThreadpoolController
controller = ThreadpoolController()


#os.environ['OPENBLAS_NUM_THREADS'] = '1' 
import math
import numpy as np
import ot
import itertools as it
from scipy.spatial.distance import *
from scipy import stats
import multiprocess
import sklearn
import shutil

from typing import Iterator, Iterable, Optional, TypeVar, Generic, Mapping


import sys



from .GW_protein import *



class Stress_Comparison:
    """
    This class streamlines computing stresses for a dataset of proteins and analyzing transport plans. 
    
    :param prot_list: A list of ``GW_protein.GW_protein`` objects
    :param RAM: Whether to store the computed transport plans in RAM versus saving to files. Default is in RAM.
    :param transport_dir: If ``RAM == False``, a filepath to save the transport plans in.
    :ivar name_list: A list of the names of the ``GW_protein``s in ``prot_list``, equivalent to ``[p.name for p in prot_list]``.   
        This serves as the keys for the dicts.
    :ivar raw_stress_dict: A dict storing all computed stresses. 
        Where ``raw_stress_dict[prot1.name][prot1.name]`` is the stress of the residues in ``prot1`` when aligning it to ``prot2``.
    :ivar dist_dict: A dict storing all computed GW or FGW distances. 
        Where ``dist_dict[prot1.name][prot1.name]`` is the GW or FGW distance between ``prot1`` and ``prot2``.
    :ivar transport_dict: A dict storing all computed transport plans; only used if ``RAM ==True``. 
        Where ``transport_dict[prot1.name][prot1.name]`` is the transport plan aligning ``prot1`` to ``prot2``.
    :ivar transport_dir: A filepath to the directory where transport plans are stored if ``RAM ==False``.
    

    """
    
    def __init__(self,
                 prot_list : list[GW_protein],
                 RAM: bool = True,
                 transport_dir: str = None):

        self.prot_list = prot_list
        self.name_list = [p.name for p in prot_list]
        self.cell_dict = {p.name: p.make_cajal_cell for p in prot_list}
        self.RAM_flag = RAM
        
        if RAM:
            self.transport_dict = {n: {n: None for n in self.name_list}  for n in self.name_list} 
        else:
            if transport_dir is None:
                raise ValueError('If RAM is not used, a directory must be provided')
                
            self.transport_dir = os.path.join(transport_dir, 'tmp_'+str(self)) 
            os.makedirs(self.transport_dir )

        
        
        
        if len(np.unique(self.name_list)) != len(self.name_list):
            raise ValueError('Names of the proteins must be unique')
        
        
        self.raw_stress_dict = {n: {n: None for n in self.name_list}  for n in self.name_list} 
        self.dist_dict = {n: {n: None for n in self.name_list}  for n in self.name_list} 
        for n in self.name_list:
            self.dist_dict[n][n] = 0
            del self.raw_stress_dict[n][n]

        self.computed_flag = False
        self.cell_dict= {p.name : p.make_cajal_cell() for p in prot_list}


    def _load_transport_plan(self,name1,name2):
        assert not self.RAM_flag
        if os.path.isfile(os.path.join(self.transport_dir, f'{name1}_{name2}.npy' )):
            T = np.load(os.path.join(self.transport_dir, f'{name1}_{name2}.npy' ) )
        elif os.path.isfile(os.path.join(self.transport_dir, f'{name2}_{name1}.npy' )):
            T = np.load(os.path.join(self.transport_dir, f'{name2}_{name1}.npy' ) ).T
        else:
            raise ValueError('transport plan not found')

        return T 


    @controller.wrap(limits=1, user_api='blas')
    def _GW_helper(self, pp):
        p1,p2 = pp 
        name1 = p1.name
        name2 = p2.name 
        cell1 = self.cell_dict[name1]
        cell2 = self.cell_dict[name2]
        c, T = GW_protein.run_GW_from_cajal(cell1, cell2, transport_plan = True)
        s1, s2 = GW_protein.GW_stress(p1,p2, T)     
        return name1, name2, c, s1, s2, T 

    @controller.wrap(limits=1, user_api='blas')
    @staticmethod
    def _GW_helper_multi(ppcc):
        pp, cc = ppcc
        p1,p2 = pp
        cell1, cell2 = cc
        name1 = p1.name
        name2 = p2.name 
        #cell1 = self.cell_dict[name1]
        #cell2 = self.cell_dict[name2]
        c, T = GW_protein.run_GW_from_cajal(cell1, cell2, transport_plan = True)
        s1, s2 = GW_protein.GW_stress(p1,p2, T)     
        return name1, name2, c, s1, s2, T 



    def GW_compute_stresses(self, processes: int = None) -> None:
        """

        This method runs all pairwise GW computations. This can be done in parallel with ``processes`` number of processes.
        
        :param processes: How many parallel processes to run, default is 1. 

        """

        
        
        
        if processes is not None and processes >1:
            with multiprocess.Pool(processes = processes)  as pool:
                cell_list = [self.cell_dict[n] for n in self.name_list]
                results = pool.imap(Stress_Comparison._GW_helper_multi, zip( it.combinations(self.prot_list,2),  it.combinations(cell_list,2)), chunksize = 20)
                for r in results:
                    name1, name2, c,s1,s2,T = r  
                    self.dist_dict[name1][name2] = c 
                    self.dist_dict[name2][name1] = c
                    self.raw_stress_dict[name1][name2] = s1
                    self.raw_stress_dict[name2][name1] = s2
                    if self.RAM_flag:
                        self.transport_dict[name1][name2] = T
                        self.transport_dict[name2][name1] = T.T
                    else:
                        np.save(file = os.path.join(self.transport_dir, f'{name1}_{name2}.npy' ), arr = T)
        else:
            results = map(self._GW_helper,   it.combinations(self.prot_list,2))
        
            for r in results:
                name1, name2, c,s1,s2,T = r  
                self.dist_dict[name1][name2] = c 
                self.dist_dict[name2][name1] = c
                self.raw_stress_dict[name1][name2] = s1
                self.raw_stress_dict[name2][name1] = s2
                if self.RAM_flag:
                    self.transport_dict[name1][name2] = T
                    self.transport_dict[name2][name1] = T.T
                else:
                    np.save(file = os.path.join(self.transport_dir, f'{name1}_{name2}.npy' ), arr = T)
        self.computed_flag = True




    @controller.wrap(limits=1, user_api='blas')
    @staticmethod
    def _FGW_lists_helper_multi(ppdda): 
        pp, dd, alpha = ppdda
        p1,p2, = pp
        data_list1, data_list2 = dd
        name1 = p1.name
        name2 = p2.name 
        n1 = len(data_list1)
        n2 = len(data_list2)

        a = np.array([np.array([x]) for x in data_list1])
        b = np.array(data_list2)
        aa = np.broadcast_to(a,(n1,n2))
        bb = np.broadcast_to(b,(n1,n2))
        M = abs(aa-bb)

        c, T = GW_protein.run_FGW_diff_mat(prot1=p1, prot2=p2, alpha = alpha, diff_mat = M , transport_plan = True)
        s1, s2 = GW_protein.FGW_stress(prot1= p1,prot2 = p2, alpha = alpha, diff_mat = M, T= T)     
        return name1, name2, c, s1, s2, T 


    def FGW_compute_stresses_data_lists(self, data_list_dict : dict, alpha: float, processes: int = None) -> None: 
        """
        This method runs all pairwise FGW computations with ``GW_protein.run_FGW_data_lists``. This can be done in parallel with ``processes`` number of processes.
        
        :param data_list_dict: A dictionary of {name: data_list} for ``GW_protein.run_FGW_data_lists``
        :param processes: How many parallel processes to run, default is 1. 
        :param alpha: The value of alpha to use for FGW.

        """

        if set(data_list_dict.keys()) != set(self.name_list):
            raise ValueError("keys of stress_list_dict must match name_list")
        for prot in self.prot_list:
            if len(prot) != len(data_list_dict[prot.name]):
                raise ValueError("lengths of data_list_dict entries must match lengths of proteins")

        data_list_list = [data_list_dict[p.name] for p in self.prot_list]
        
        if processes is not None and processes >1:
            with multiprocess.Pool(processes = processes)  as pool:
                results = pool.imap(Stress_Comparison._FGW_lists_helper_multi, zip( it.combinations(self.prot_list,2), it.combinations(data_list_list,2), it.repeat(alpha)), chunksize = 20)
                for r in results:
                    name1, name2, c,s1,s2,T = r  
                    self.dist_dict[name1][name2] = c 
                    self.dist_dict[name2][name1] = c
                    self.raw_stress_dict[name1][name2] = s1
                    self.raw_stress_dict[name2][name1] = s2
                    if self.RAM_flag:
                        self.transport_dict[name1][name2] = T
                        self.transport_dict[name2][name1] = T.T
                    else:
                        np.save(file = os.path.join(self.transport_dir, f'{name1}_{name2}.npy' ), arr = T)
        else:
            results = map(Stress_Comparison._FGW_lists_helper_multi, zip( it.combinations(self.prot_list,2), it.combinations(data_list_list,2), it.repeat(alpha)))
        
            for r in results:
                name1, name2, c,s1,s2,T = r  
                self.dist_dict[name1][name2] = c 
                self.dist_dict[name2][name1] = c
                self.raw_stress_dict[name1][name2] = s1
                self.raw_stress_dict[name2][name1] = s2
                if self.RAM_flag:
                    self.transport_dict[name1][name2] = T
                    self.transport_dict[name2][name1] = T.T
                else:
                    np.save(file = os.path.join(self.transport_dir, f'{name1}_{name2}.npy' ), arr = T)
        self.computed_flag = True


    @controller.wrap(limits=1, user_api='blas')
    @staticmethod
    def _FGW_dict_helper_multi(ppda): 
        pp, ddict, alpha = ppda
        p1,p2, = pp
        name1 = p1.name
        name2 = p2.name 
        n1 = len(p1)
        n2 = len(p2)
        M = np.zeros((n1,n2))
        for i in range(n1):
            for j in range(n2):
                M[i,j] = ddict[p1.seq[i]][p2.seq[j]]

        c, T = GW_protein.run_FGW_diff_mat(prot1=p1, prot2=p2, alpha = alpha, diff_mat = M, transport_plan = True)
        s1, s2 = GW_protein.FGW_stress(prot1= p1,prot2 = p2, alpha = alpha, diff_mat = M, T= T)     
        return name1, name2, c, s1, s2, T 
        
    def FGW_compute_stresses_dict(self, diff_dict : dict, alpha: float, processes: int = None) -> None: 
        """

        This method runs all pairwise FGW computations with ``GW_protein.run_FGW_dict``. This can be done in parallel with ``processes`` number of processes.
        
        :param data_list_dict: A dictionary for ``GW_protein.run_FGW_dict``
        :param processes: How many parallel processes to run, default is 1. 
        :param alpha: The value of alpha to use for FGW.

        """


        if processes is not None and processes >1:
            with multiprocess.Pool(processes = processes)  as pool:
                results = pool.imap(Stress_Comparison._FGW_dict_helper_multi, zip( it.combinations(self.prot_list,2), it.repeat(diff_dict), it.repeat(alpha)), chunksize = 20)
                for r in results:
                    name1, name2, c,s1,s2,T = r  
                    self.dist_dict[name1][name2] = c 
                    self.dist_dict[name2][name1] = c
                    self.raw_stress_dict[name1][name2] = s1
                    self.raw_stress_dict[name2][name1] = s2
                    if self.RAM_flag:
                        self.transport_dict[name1][name2] = T
                        self.transport_dict[name2][name1] = T.T
                    else:
                        np.save(file = os.path.join(self.transport_dir, f'{name1}_{name2}.npy' ), arr = T)
        else:
            results = map(Stress_Comparison._FGW_lists_helper_multi, zip( it.combinations(self.prot_list,2), it.repeat(diff_dict), it.repeat(alpha)))
        
            for r in results:
                name1, name2, c,s1,s2,T = r  
                self.dist_dict[name1][name2] = c 
                self.dist_dict[name2][name1] = c
                self.raw_stress_dict[name1][name2] = s1
                self.raw_stress_dict[name2][name1] = s2
                if self.RAM_flag:
                    self.transport_dict[name1][name2] = T
                    self.transport_dict[name2][name1] = T.T
                else:
                    np.save(file = os.path.join(self.transport_dir, f'{name1}_{name2}.npy' ), arr = T)
        self.computed_flag = True


    
    def __del__(self):
        #print('__del__ called', self) #debugging
        if not self.RAM_flag:
            shutil.rmtree(self.transport_dir)




    def raw_transferred_stresses(self, stress_dict: dict[str,np.array]) ->dict[str,dict[str,np.array]]:
        """

        This method computes all of the transferred stresses. 
        
        :param stress_dict: A dictionary of the form ``{name: np.array}`` with keys the protein names in ``self.name_list``, where the arrays represent the stresses of each residue.
        :return dict: A dictionary of all transferred stresses where
                Where ``raw_transferred_stress[prot1.name][prot1.name]`` is the stress of the residues in ``prot1`` based on the transferred stress from ``prot2``.


        """
        assert set(stress_dict.keys()) == set(self.name_list)
        assert self.computed_flag
        for i in range(len(self.prot_list)):
            assert len(self.prot_list[i]) == len(stress_dict[self.name_list[i]])
        raw_transferred_stresses = {n: {n: None for n in self.name_list}  for n in self.name_list} 
        for n in self.name_list:
            del raw_transferred_stresses[n][n]
        
        
        for nn in it.combinations(self.name_list,2):
            name1, name2 = nn
            if self.RAM_flag:
                T = self.transport_dict[name1][name2]
            else:
                T = self._load_transport_plan(name1,name2)
        
            raw_transferred_stresses[name1][name2] = T @ stress_dict[name2]
            raw_transferred_stresses[name2][name1] = T.T@ stress_dict[name1]
        
        
        return raw_transferred_stresses


    def get_GW_dmat(self) -> np.array:
        """

        This method returns the GW or FGW distances of ``self`` in the form of a distance matrix. 
        The indexing is that of ``self.prot_list``
        
        :return: np.array
        
        """

        assert self.computed_flag

        n = len(self.name_list)
        dmat = np.zeros((n,n))
        for i in range(n):
            for j in range(i):
                name1 = self.name_list[i]
                name2 = self.name_list[j]
                assert self.dist_dict[name1][name2] == self.dist_dict[name2][name1]
                if self.dist_dict[name1][name2] is None:
                    raise RuntimeError('Stresses and distances must be computed with \
                        GW_compute_stresses before get_GW_dmat can be run')
                else:
                    dmat[i,j] = self.dist_dict[name1][name2]
                    dmat[j,i] = self.dist_dict[name1][name2]
        return dmat




def normalize_stress_dict(raw_dict: dict[str,dict[str,np.array]], code: tuple[float,float,float,float] =(1, 0, 0, 0, 0)) -> dict[str,np.array]: 
    """
    This method takes in a dictionary of raw stress and outputs a dictionary of weighted averages.

    :param raw_dict: A dictionary of raw stresses of the format ``raw_dict[name1][name2] == stress``.
    :param code: This is a tuple of exponents to be used for weighting.
        ``code[0]`` is the exponent for each stress value, 
        ``code[1]`` is the exponent for each stress value to be summed, 
        ``code[2]`` is the exponent of the total stress in a row,
        ``code[3]`` is the exponent of the number of residues in the protein,
        ``code[4]`` is the exponent for the number of other proteins.
        The default is (1,0,0,0,0) which corresponds to the simple sum. (1,0,0,0,-1) is the mean. 
    :return dict: A dictionary of stresses of the format ``stress_dict[name]== stress``
    
    """
    a, b, c, d, e = code
    norm_stresses_dict = {}
    for k in raw_dict.keys():
        mat = np.stack(list(raw_dict[k].values())) 
        out = np.sum(
            mat**a * (np.sum(mat**b, axis=1) ** c)[:, np.newaxis] * mat.shape[1] ** d,
            axis=0,
        )* mat.shape[0] ** e
        norm_stresses_dict[k] = out
    return norm_stresses_dict


def get_percentile_of_dict(stress_dict: dict[str,np.array]) ->dict[str,np.array]:
    """
    This method replaces each stress array with the percentile values of the stresses for that protein.

    :param stress_dict: A dict of the stress levels for each protein.
    :return: A dict of the percentiles of the stress levels at each residue for each protein.

    """
    return {
        k: scipy.stats.percentileofscore(stress_dict[k], stress_dict[k])
        for k in stress_dict.keys()
    }


def get_AP_scores(stress_dict: dict[str,np.array], true_region_dict : dict[str,list[bool]], upper = False ) ->dict[str, float]:
    """
    This method takes a stress dict and calculates the average precision for each protein
     of using the stresses to predict user-inputted regions in the proteins

    :param stress_dict: A dict of the stress levels for each protein.
    :param true_region_dict: A dict of the true regions to be predicted.
    :param upper: Whether to predict the regions based on high stress (``True``) or low stress (``False``)
    :return: A dict of the average precision scores for each protein.

    """

    if upper:
            AP_dict = {
            name: sklearn.metrics.average_precision_score(
                y_true=true_region_dict[name], y_score=[ s for s in stress_dict[name]]
            )
            for name in stress_dict.keys()
        }
    else:
        AP_dict = {
            name: sklearn.metrics.average_precision_score(
                y_true=true_region_dict[name], y_score=[1 - s for s in stress_dict[name]]
            )
            for name in stress_dict.keys()
        }
    return AP_dict


# def mean_AP_scores(stress_dict: dict[str,np.array], true_region_dict : dict[str,list[bool]],upper = False ):
#     if upper:
#         return np.mean(
#             [
#                 metrics.average_precision_score(
#                     y_true=true_region_dict[name],
#                     y_score=[ s for s in stress_dict[name]],
#                 )
#                 for name in stress_dict.keys()
#             ]
#         )
        
#     else:
#         return np.mean(
#             [
#                 metrics.average_precision_score(
#                     y_true=true_region_dict[name],
#                     y_score=[1 - s for s in stress_dict[name]],
#                 )
#                 for name in stress_dict.keys()
#             ]
#         )


# def single_threshold_AP_score(stress_dict: dict[str,np.array], true_region_dict : dict[str,list[bool]],  upper = False ):
#     full_stresses = []
#     full_true_regions = []
#     for name in list(stress_dict.keys()):
#         full_stresses += list(stress_dict[name])
#         full_true_regions += list(true_region_dict[name])

#     if upper:
#         return metrics.average_precision_score(
#         y_true=np.array(full_true_regions),
#         y_score=np.array([ s for s in full_stresses]),
#         )
#     else:
#         return metrics.average_precision_score(
#             y_true=np.array(full_true_regions),
#             y_score=np.array([1 - s for s in full_stresses]),
#         )


# def avgd_single_threshold_AP_score(stress_dict: dict[str,np.array], true_region_dict : dict[str,list[bool]],  upper = False ):
#     avgd_full_stresses = []
#     full_true_regions = []
#     for name in list(stress_dict.keys()):
#         avgd_full_stresses += list(
#             stress_dict[name] / np.mean(stress_dict[name])
#         )
#         full_true_regions += list(true_regions_dict[name])

#     if upper:
#         return metrics.average_precision_score(
#             y_true=np.array(full_true_regions),
#             y_score=np.array([s for s in avgd_full_stresses]),
#         )
#     else:
#         return metrics.average_precision_score(
#             y_true=np.array(full_true_regions),
#             y_score=np.array([1 - s for s in avgd_full_stresses]),
#         )


# def z_single_threshold_AP_score(stress_dict: dict[str,np.array], true_region_dict : dict[str,list[bool]],  upper = False) ->float:
#     z_full_stresses = []
#     full_true_regions = []
#     for name in list(stress_dict.keys()):
#         z_full_stresses += list(scipy.stats.zscore(stress_dict[name]))
#         full_true_regions += list(true_region_dict[name])

#     if upper:
#         return metrics.average_precision_score(
#             y_true=np.array(full_true_regions),
#             y_score=np.array([s for s in z_full_stresses]),
#         )
#     else:
#         return metrics.average_precision_score(
#             y_true=np.array(full_true_regions),
#             y_score=np.array([1 - s for s in z_full_stresses]),
#         )



