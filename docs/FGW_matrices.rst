
FGW_matrices
==========================

This module contains methods to generate difference dicts and data lists for the fused Gromov-Wasserstein methods in  ``GWProt.GW_protein.run_FGW_dict`` and ``GWProt.GW_protein.run_FGW_data_lists``.

.. module:: GWProt.FGW_matrices


----------------------




.. autofunction:: GWProt.FGW_matrices.get_hydrophobicity_dict

.. autofunction:: GWProt.FGW_matrices.get_hydrophobicity_list

These methods give difference matrices using the hydrophobicity values in `Eisenberg, Schwarz, Komaromy, and Wall <https://pubmed.ncbi.nlm.nih.gov/6502707/>`_.

----------------------------

.. autofunction:: GWProt.FGW_matrices.get_BLOSUM_dict


This method gives difference matrices based on the `BLOSUM matrices <https://en.wikipedia.org/wiki/BLOSUM>`_. For a pair of amino acids it first takes :math:`e^{-b}`, where :math:`b` is the BLOSUM entry, then normalizes to make the distances satisfy the triangle inequality.

----------------------------



.. autofunction:: GWProt.FGW_matrices.get_Grantham_dict

This method gives difference matrices given by the `Grantham <https://www.jstor.org/stable/1739007?seq=1>`_ differnce scores.

--------------------



.. autofunction:: GWProt.FGW_matrices.get_pI_dict

.. autofunction:: GWProt.FGW_matrices.get_pI_list

These methods give difference matrices using the Solomon isoelectric point values. It is preferable to use the ``GWProt.GW_pI`` module instead as it has greater functionality.