.. GWProt documentation master file


Welcome to GWProt's documentation
===================================

**GWProt** is a computational framework that leverages recent advances in metric geometry, specifically the use of Gromov-Wasserstein (GW) 
correspondences, for protein structure alignment. GW correspondences find an optimal structural alignment between pairs of proteins by 
minimizing the distortion of intra-molecular distances among paired residues. GWProt enables the incorporation of biochemical information 
into structural comparisons and introduces the concept of local geometric distortion, a measure that captures fine-scale conformational 
differences. Using this framework, we can identify conformational switches in individual proteins, detect functional domains shared among 
evolutionarily distant proteins, reveal topological rearrangements in homologous folds, or uncover recurrent short structural motifs 
underlying functional domains. GWProt has been developed by the `Cámara Lab <https://camara-lab.org/>`_. Installation instructions can be 
found at the `GWProt GitHub repository <https://github.com/CamaraLab/GWProt>`_.


.. toctree::
   :maxdepth: 1
   :caption: Background

   gromov_wasserstein
   fused_gromov_wasserstein
   stress
   distortion_scaling


.. toctree::
   :maxdepth: 1
   :caption: Examples

   Examples/Example_1_Comparing_KRAS_Proteins
   Examples/Example_2_Permutations
   Examples/Example_3_Clustering_RdRps

.. toctree::
   :maxdepth: 2
   :caption: API

   GW_protein
   GW_protein_pI
   FGW_matrices
   stress_comparison
   switch_probabilities
   pymol_protein_viewer
   weighted_alignment
   




Indices
==================

* :ref:`genindex`
* :ref:`search`