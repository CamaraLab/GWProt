The GW_protein Class
==========================

This class has the core functionalities of GWProt

.. module:: GWProt.GW_protein
.. autoclass:: GWProt.GW_protein.GW_protein


A `GW_protein` object contains all the the data used to compute the GW distance.


------------------


# Basic Methods
We have basic ways to create and compare `GW_protein` objects.

The usual way to make a `GW_protein` object is by loading it from a ``.pdb`` (Protein Data Bank) file. 
.. autofunction:: GWProt.GW_protein.GW_protein.make_protein_from_pdb

If data is missing in the form of missing residues or missing alpha-Carbons it will be skipped. 

.. autofunction:: GWProt.GW_protein.GW_protein.validate
.. autofunction:: GWProt.GW_protein.GW_protein.__eq__
.. autofunction:: GWProt.GW_protein.GW_protein.__len__

------------------

# Intra-Protein Distance Matrix Manipulation
Next we have methods to manipulate the intra-protein distance matrix.
`distortion_scaling`_

.. autofunction:: GWProt.GW_protein.GW_protein.scale_ipdm
.. autofunction:: GWProt.GW_protein.GW_protein.reset_ipdm

------------------
# Downsampling
Then we have two methods for downsampling. 
Downsampling reduces the number of residues used so has the effect of speeding up computations, but can reduce accuracy.

.. autofunction:: GWProt.GW_protein.GW_protein.downsample_by_indices
.. autofunction:: GWProt.GW_protein.GW_protein.downsample_n



-------------------
# Computing GW

 `gromov_wasserstein`_

The methods for computing the Gromov-Wasserstein distance use the [CAJAL library](https://github.com/CamaraLab/CAJAL/tree/main), also created by the CamaraLab, for efficient computation. 


.. autofunction:: GWProt.GW_protein.GW_protein.make_GW_cell
.. autofunction:: GWProt.GW_protein.GW_protein.run_GW_from_cells
.. autofunction:: GWProt.GW_protein.GW_protein.run_GW


As this uses CAJAL, there is the ability to use other functionalities from CAJAL.

--------------------
# Computing FGW


`fused_gromov_wasserstein`_


Multiple ways of inputting the feature space data are included.

.. autofunction:: GWProt.GW_protein.GW_protein.run_FGW_diff_mat

The first is the most general as it can use any user-inputted feature difference matrix. However a new difference matrix must be used for every pair of proteins. 

.. autofunction:: GWProt.GW_protein.GW_protein.run_FGW_data_lists

The second uses a linear feature space. This is suitable for scalar features including isoelectric point, solvent-accessible surface area, charge, and hydrophobicity. 


.. autofunction:: GWProt.GW_protein.GW_protein.run_FGW_dict

The third uses 
The input ``dict`` 




It is not recommended to use these on downsampled proteins, as the data is lost from the excluded residues. 

As CAJAL does not run fused GW, these computation are done with the Python ``ot`` library. 


--------------------

# Computing Stress 

`stress`_


.. autofunction:: GWProt.GW_protein.GW_protein.GW_stress
.. autofunction:: GWProt.GW_protein.GW_protein.FGW_stress

WARNING - ``np.sum(stress1) != c``, where ``c`` is the GW cost; rather ``math.sqrt(np.sum(stress1))/2 == c``; and similarly for ``stress2``.


------------------

# Miscellaneous Methods

.. autofunction:: GWProt.GW_protein.GW_protein.get_eccentricity
.. autofunction:: GWProt.GW_protein.GW_protein.run_ssearch_indices




