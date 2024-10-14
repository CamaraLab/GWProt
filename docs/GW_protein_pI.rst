The GW_protein_pI Class
==========================

.. module:: GWProt.GW_protein_pI

.. autoclass:: GWProt.GW_protein_pI.GW_protein_pI


Unlike ``GW_protein`` , this stores the isoelectric point values with each protein [1]. Most of the methods are the same as those of ``GW_protein``. A key difference is that when downsampling, we can combine the isoelectric points of adjacent residues which are grouped together, giving an estimated isoelectric point of the segments. Thus unlike the FGW methods in ``GW_protein`` , we can run FGW on downsampled proteins in a meaningful way. This is done using an algorithm based on the one in the `Sequence Manipulation Suite <https://github.com/paulstothard/sequence_manipulation_suite>`_ based on the Henderson-Hasselbach equation.

For uniform downsampling, this is done automatically in the ``downsample_n`` method.

.. autofunction:: GWProt.GW_protein_pI.GW_protein_pI.downsample_n



For downsampling to specified indices, we first need to ' smooth out ' the isoelectic values, so that when we downsample, we will be accounting for the isoelectric points of nearby discarded residues.
 
.. autofunction:: GWProt.GW_protein_pI.GW_protein_pI.convolve_pIs


---------------------

Then we have the same methods as in ``GW_protein`` though adapted to use isoelectric points.


.. autofunction:: GWProt.GW_protein_pI.GW_protein_pI.validate


.. autofunction:: GWProt.GW_protein_pI.GW_protein_pI.__eq__



.. autofunction:: GWProt.GW_protein_pI.GW_protein_pI.downsample_by_indices


.. autofunction:: GWProt.GW_protein_pI.GW_protein_pI.run_FGW


.. autofunction:: GWProt.GW_protein_pI.GW_protein_pI.FGW_stress


.. autofunction:: GWProt.GW_protein_pI.GW_protein_pI.run_FGW_seq_aln


.. [1] We note that these are rather naive estimates of the isoelectric points. More sophisticated ones can be computed by other software packages using the 3-dimensional structure of a protein. This could then be used with ``GW_protein.GW_protein.run_FGW_data_lists`` .

