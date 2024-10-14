GW_protein
==========================

This class has the core functionalities of GWProt. A ``GW_protein`` object contains all the the data used to compute the GW distance.


.. module:: GWProt.GW_protein
.. autoclass:: GWProt.GW_protein.GW_protein



**Basic Methods**
---------------------



We have basic ways to create and compare ``GW_protein`` objects.

The usual way to make a ``GW_protein`` object is by loading it from a ``.pdb`` (Protein Data Bank) file. 

.. autofunction:: GWProt.GW_protein.GW_protein.make_protein_from_pdb

If data is missing in the form of missing residues or missing alpha-Carbons it will be skipped. Note that all indices within a ``GW_protein`` object are based on those loaded, which may not agree with the indices in the pdb file. 

.. autofunction:: GWProt.GW_protein.GW_protein.validate

.. autofunction:: GWProt.GW_protein.GW_protein.__eq__

.. autofunction:: GWProt.GW_protein.GW_protein.__len__


**Intra-Protein Distance Matrix Scaling**
--------------------------------------------------


Next we have methods to manipulate the intra-protein distance matrix for distortion scaling.

.. autofunction:: GWProt.GW_protein.GW_protein.scale_ipdm

.. autofunction:: GWProt.GW_protein.GW_protein.reset_ipdm


**Downsampling**
-------------------


Then we have two methods for downsampling. 
Downsampling reduces the number of residues used so has the effect of speeding up computations, but can reduce accuracy.

.. autofunction:: GWProt.GW_protein.GW_protein.downsample_by_indices

.. autofunction:: GWProt.GW_protein.GW_protein.downsample_n



**Computing GW**
-------------------



The methods for computing the Gromov-Wasserstein distance use the `CAJAL library <https://github.com/CamaraLab/CAJAL/tree/main>`_ , also created by the CámaraLab, for efficient computation. 


.. autofunction:: GWProt.GW_protein.GW_protein.make_GW_cell

.. autofunction:: GWProt.GW_protein.GW_protein.run_GW_from_cells

.. autofunction:: GWProt.GW_protein.GW_protein.run_GW


As this uses CAJAL, there is the ability to use other functionalities from CAJAL.

**Computing FGW**
--------------------





Multiple ways of inputting the feature space data are included.

.. autofunction:: GWProt.GW_protein.GW_protein.run_FGW_diff_mat

The first is the most general as it can use any user-inputted feature difference matrix. However a new difference matrix must be used for every pair of proteins. 

.. autofunction:: GWProt.GW_protein.GW_protein.run_FGW_data_lists

The second uses a linear feature space. This is suitable for scalar features including isoelectric point, solvent-accessible surface area, charge, and hydrophobicity. 


.. autofunction:: GWProt.GW_protein.GW_protein.run_FGW_dict

The third uses a dictionary giving difference values between different types of amino acids.




It is not recommended to use these on downsampled proteins, as the data is lost from the excluded residues. 

As CAJAL does not run FGW, these computation are done with the Python ``ot`` library. 


**Computing Stress**
--------------------




.. autofunction:: GWProt.GW_protein.GW_protein.GW_stress

.. autofunction:: GWProt.GW_protein.GW_protein.FGW_stress

WARNING - ``np.sum(stress1) != c``, where ``c`` is the GW cost; rather ``math.sqrt(np.sum(stress1))/2 == c``; and similarly for ``stress2`` , and for FGW.


**Miscellaneous Methods**
---------------------------------

.. autofunction:: GWProt.GW_protein.GW_protein.run_ssearch_indices

Explicity this runs the command

``$ ssearch36 -s BP62 -p -T 1 -b 1 -f 0 -g 0 -z -1 -m 9C``

for the Smith-Waterman algorithm in the `Fasta Package <https://github.com/wrpearson/fasta36>`_.

.. autofunction:: GWProt.GW_protein.GW_protein.run_GW_seq_aln




.. autofunction:: GWProt.GW_protein.GW_protein.get_eccentricity

 Eccentricity is defined in `Memolis's paper <https://www.math.ucdavis.edu/~saito/data/acha.read.w12/memoli-gromov-dist.pdf>`_, Definition 5.3. Intuitively it quantifies how far each residue is from the rest of the residues in a protein. Within a give protein, residues with higher eccentricity often have higher stress when aligned to other proteins, so this could be used for normalization.





