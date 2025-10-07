
FGW_matrices
==========================


This module provides methods to generate difference dictionaries and data lists for the fused Gromov-Wasserstein (FGW) methods in ``GWProt.GW_protein.run_FGW_dict`` and ``GWProt.GW_protein.run_FGW_data_lists``. These difference matrices are used to incorporate biochemical features into FGW alignments and to compute local geometric distortion (LGD) at the residue level.

.. module:: GWProt.FGW_matrices


----------------------





**Hydrophobicity-Based Differences**
-----------------------------------
.. autofunction:: GWProt.FGW_matrices.get_hydrophobicity_dict
.. autofunction:: GWProt.FGW_matrices.get_hydrophobicity_list

These methods generate difference matrices using hydrophobicity values from `Eisenberg, Schwarz, Komaromy, and Wall <https://pubmed.ncbi.nlm.nih.gov/6502707/>`_. These can be used as the feature space for FGW, allowing LGD to reflect hydrophobicity differences between residues.

----------------------------


**BLOSUM-Based Differences**
---------------------------
.. autofunction:: GWProt.FGW_matrices.get_BLOSUM_dict

This method generates difference matrices based on the `BLOSUM matrices <https://en.wikipedia.org/wiki/BLOSUM>`_. For a pair of amino acids, it computes :math:`e^{-b}` (where :math:`b` is the BLOSUM entry), then normalizes the result to ensure the distances satisfy the triangle inequality. These matrices can be used in FGW to compute LGD based on sequence similarity.

----------------------------




**Grantham-Based Differences**
-----------------------------
.. autofunction:: GWProt.FGW_matrices.get_Grantham_dict

This method generates difference matrices using the `Grantham <https://www.jstor.org/stable/1739007?seq=1>`_ difference scores, which reflect physicochemical differences between amino acids. These can be used in FGW to compute LGD based on these properties.

--------------------




**Isoelectric Point Differences**
---------------------------------
.. autofunction:: GWProt.FGW_matrices.get_pI_dict
.. autofunction:: GWProt.FGW_matrices.get_pI_list

These methods generate difference matrices using Solomon isoelectric point values. For more advanced handling of isoelectric points, use the ``GWProt.GW_pI`` module, which offers greater functionality for FGW and LGD computations.