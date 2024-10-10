
The FGW_matrices
==========================

This module contains methods to generate difference matrices for the fused Gromov-Wasserstein methods in  ``GWProt.GW_protein.run_FGW_dict`` and ``GWProt.GW_protein.run_FGW_data_lists``.


These methods give difference matrices using the hydrophobicity values in [Eisenberg, Schwarz, Komaromy, and Wall](https://pubmed.ncbi.nlm.nih.gov/6502707/).
.. autofunction:: GWProt.FGW_matrices.get_hydrophobicity_dict
.. autofunction:: GWProt.FGW_matrices.get_hydrophobicity_list


This method gives difference matrices based on the [BLOSUM matrices](https://en.wikipedia.org/wiki/BLOSUM). For a pair of amino acids it first takes :math: `e^{-b}`, where :math:`b` is the BLOSUM entry, then normalizes to make the distances satisfy the triangle inequality.
.. autofunction:: GWProt.FGW_matrices.get_BLOSUM_dict



This method gives difference matrices given by the [Grantham](https://www.jstor.org/stable/1739007?seq=1) differnce scores.
.. autofunction:: GWProt.FGW_matrices.get_Grantham_dict


These methods give difference matrices using the Soloman isoelectric point values. It is preferable to use the ``GWProt.GW_pI`` module instead as it has greater functionality.
.. autofunction:: GWProt.FGW_matrices.get_PI_dict
.. autofunction:: GWProt.FGW_matrices.get_PI_list
