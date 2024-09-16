import time
import re
import os

os.environ['OPENBLAS_NUM_THREADS'] = '1' 
import math
import numpy as np
import random
import ot
import statistics
import itertools as it
from scipy.spatial.distance import *
import multiprocessing
import multiprocess
import pandas as pd
import shutil
import cajal


import Bio.PDB
from typing import Iterator, Iterable, Optional, TypeVar, Generic

# imports used in compute_gw_copy.ipynb
from Bio.SVDSuperimposer import SVDSuperimposer


import sys



import GW_scripts
import read_pdb
import FGW_protein





# This module is for doing all vs all stress comparisons in a dataset

class Stress_Comparison:
    
    def __init__(self, prot_list, RAM = True, transport_dir = None):
        self.prot_list = prot_list
        self.name_list = [p.name for p in prot_list]
        self.cell_dict = {p.name: p.make_GW_cell for p in prot_list}
        self.RAM_flag = RAM
        
        if RAM:
            self.transport_dict = {n: {n: None for n in self.name_list}  for n in self.name_list} 
        else:
            if transport_dir is None:
                raise ValueError('If RAM is not used, a directory must be provided')
            if not os.path.isdir(transport_dir):
                raise ValueError(f'The directory {transport_dir} is not found')
        
        
        
        if len(np.unique(self.name_list)) != len(self.name_list):
            raise ValueError('Names of the proteins must be unique')
        self.transport_dir = transport_dir
        
        
        self.raw_stress_dict = {n: {n: None for n in self.name_list}  for n in self.name_list} 
        self.dist_dict = {n: {n: None for n in self.name_list}  for n in self.name_list} 
        for n in self.name_list:
            self.dist_dict[n][n] = 0
        
        self.cell_dict= {p.name : p.make_GW_cell() for p in prot_list}
            #makes them with uniform distribution


    def _load_transport_plan(name1,name2):
        assert not self.RAM_flag
        if os.path.isfile(os.path.join(transport_dir, f'{name1}_{name2}.npy' )):
            T = np.load(os.path.join(transport_dir, f'{name1}_{name2}.npy' ) )
        elif os.path.isfile(os.path.join(transport_dir, f'{name2}_{name1}.npy' )):
            T = np.load(os.path.join(transport_dir, f'{name2}_{name1}.npy' ) ).T
        else:
            raise ValueError('transport plan not found')

        return T 



    def _GW_helper(self, pp):
        p1,p2 = pp 
        name1 = p1.name
        name2 = p2.name 
        cell1 = self.cell_dict[name1]
        cell2 = self.cell_dict[name2]
        c, T = FGW_protein.FGW_protein.run_GW_from_cells(cell1, cell2, transport_plan = True)
        s1, s2 = FGW_protein.FGW_protein.GW_stress(p1,p2, T)     
        print(name1, name2, c)
        return name1, name2, c, s1, s2, T 

    @staticmethod
    def _GW_helper_multi(ppcc):
        pp, cc = ppcc
        p1,p2 = pp
        cell1, cell2 = cc
        name1 = p1.name
        name2 = p2.name 
        #cell1 = self.cell_dict[name1]
        #cell2 = self.cell_dict[name2]
        c, T = FGW_protein.FGW_protein.run_GW_from_cells(cell1, cell2, transport_plan = True)
        s1, s2 = FGW_protein.FGW_protein.GW_stress(p1,p2, T)     
        print(name1, name2, c)
        return name1, name2, c, s1, s2, T 


    def GW_compute_stresses(self, processes = None):
        #core None, 0, or 1 makes it not multiprocess
        #otherwise it gives the number of cores
        
        #todo - add progress bar
        
        
        if processes is not None and processes >1:
            with multiprocessing.Pool(processes = processes)  as pool:
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

    
    def __del__(self):
        print('__del__ called', self) #debugging
        if not self.RAM_flag:
            shutil.rmtree(self.transport_dir)


    def _stress_transfer_helper(self, stress_dict, pp):
        p1, p2 = pp
        name1 = p1.name
        name2 = p2.name
        if self.RAM_flag:
            T = self.transport_dict[name1][name2]
        else:
            T = _load_transport_plan(name1,name2)

        s1 = T @ stress_dict[name2]
        s2 = T.T@ stress_dict[name1]
        return name1, name2, s1,s2


    def raw_transferred_stresses(self, stress_dict, cores = None):
        assert set(stress_dict.keys()) == set(self.name_list)
        raw_transferred_stresses = {n: {n: None for n in name_list}  for n in name_list} 

         

        if cores is not None and cores >1:
            with multiprocessing.Pool(cores) as pool:
                results = pool.imap(helper, zip( it.repeat(stress_dict), it.combinations(prot_list)), chunksize = 20)
        else:
            results = map(helper, zip( it.repeat(stress_dict), it.combinations(prot_list)))

        for r in results:
            n1,n2, s1,s2 = r
            raw_transferred_stresses[n1][n2] = s1 
            raw_transferred_stresses[n2][n1] = s2
        return raw_transferred_stresses


    def get_GW_dmat(self):
        #returns a np.array of the GW distances
        #ordering is that of the prot_list
        n = len(name_list)
        dmat = np.zeros((n,n))
        for i in range(n):
            for j in range(i):
                name1 = name_list[i]
                name2 = name_list[j]
                assert dist_dict[name1][name2] == dist_dict[name2][name1]
                if dist_dict[name1][name2] is None:
                    raise RuntimeError('Stresses and distances must be computed with \
                        GW_compute_stresses before get_GW_dmat can be run')
                else:
                    dmat[i,j] = dist_dict[name1][name2]
                    dmat[j,i] = dist_dict[name1][name2]
        return dmat















def normalize_stress_dict(raw_dict, code=(1, 0, 0, 0, 0)):
    a, b, c, d, ecc = code
    norm_stresses_dict = {}
    for k in raw_dict.keys():
        mat = np.array(raw_dict[k])
        out = np.sum(
            mat**a * (np.sum(mat**b, axis=1) ** c)[:, np.newaxis] * mat.shape[1] ** d,
            axis=0,
        )
        norm_stresses_dict[k] = out
    return norm_stresses_dict


def get_percentile_of_dict(stress_dict):
    return {
        k: scipy.stats.percentileofscore(stress_dict[k], stress_dict[k])
        for k in stress_dict.keys()
    }


def get_AP_scores(stress_dict, true_region_dict ):
    AP_dict = {
        name: sklearn.metrics.average_precision_score(
            y_true=true_region_dict[name], y_score=[1 - s for s in stress_dict[name]]
        )
        for name in stress_dict.keys()
    }
    return AP_dict


def mean_AP_scores(stress_dict, true_region_dict ):
    return np.mean(
        [
            metrics.average_precision_score(
                y_true=true_region_dict[name],
                y_score=[1 - s for s in stress_dict[name]],
            )
            for name in stress_dict.keys()
        ]
    )


def single_threshold_AP_score(normalized_stress_dict, true_region_dict ):
    full_stresses = []
    full_true_regions = []
    for name in list(stress_dict.keys()):
        full_stresses += list(normalized_stress_dict[name])
        full_true_regions += list(true_region_dict[name])
    return metrics.average_precision_score(
        y_true=np.array(full_true_regions),
        y_score=np.array([1 - s for s in full_stresses]),
    )


def avgd_single_threshold_AP_score(normalized_stress_dict, true_region_dict ):
    avgd_full_stresses = []
    full_true_regions = []
    for name in list(stress_dict.keys()):
        avgd_full_stresses += list(
            normalized_stress_dict[name] / np.mean(normalized_stress_dict[name])
        )
        full_true_regions += list(true_regions_dict[name])
    return metrics.average_precision_score(
        y_true=np.array(full_true_regions),
        y_score=np.array([1 - s for s in avgd_full_stresses]),
    )


def z_single_threshold_AP_score(normalized_stress_dict, true_region_dict ):
    z_full_stresses = []
    full_true_regions = []
    for name in list(stress_dict.keys()):
        z_full_stresses += list(scipy.stats.zscore(normalized_stress_dict[name]))
        full_true_regions += list(true_region_dict[name])
    return metrics.average_precision_score(
        y_true=np.array(full_true_regions),
        y_score=np.array([1 - s for s in z_full_stresses]),
    )





def get_single_prec_rec(stress_dict,  true_region_dict,threshold=0.5 ):
    full_stresses = []
    full_true_regions = []
    for name in list(stress_dict.keys()):
        full_stresses += list(stress_dict[name])
        full_true_regions += list(true_region_dict[name])

    TP = 0
    FP = 0
    FN = 0
    for j in range(len(full_true_regions)):
        if full_true_regions[j] and full_stresses[j] <= threshold:  # TP
            TP += 1
            # if full_stresses[j] == threshold:
            #     print(' +, == threshold')
        elif full_true_regions[j] and full_stresses[j] > threshold:  # FN
            FN += 1
        elif not full_true_regions[j] and full_stresses[j] < threshold:  # FP
            FP += 1
        elif not full_true_regions[j] and full_stresses[j] == threshold:  # FP
            FP += 1
            # print(' - ,== threshold')
    return (
        TP / (TP + FP),
        TP / (TP + FN) if (TP + FN) * (TP + FP) > 0 else -1,
        TP,
        FN,
        FP,
    )

