.. GWProt documentation master file


Welcome to GWProt's documentation!
===================================

**GWProt**  is a Python library for applying Gromov-Wasserstein to protein morphology developed by the `Cámara Lab <https://camara-lab.org/>`_



.. note::

   This project is under active development. Documentation and code will be updated continuously


With hundreds of thousands of protein structures experimentally determined and hundreds of millions predicted, computational methods for comparing protein structures are becoming key to determining evolutionary connections and functional similarities. Most existing methods rely on linear sequence alignment and rigid structural alignment, which limits their utility in studying circularly permuted and other multiply rearranged proteins, especially those with low sequence homology. 
We introduce GWProt: a novel method of protein structure alignment based on concepts from metric geometry and optimal transport theory. GWProt utilizes the Gromov-Wasserstein (GW) distance to find an optimal structural alignment between pairs of proteins that minimizes the distortion of intra-molecular distances among paired residues. In addition, the flexibility of GW distances presents some unique features, including the possibility of incorporating biochemical information such as isoelectric points into the alignment, as well as user-inputted data.


**Installation:**
GWProt can be installed by running

``pip install GWProt@git+https://github.com/CamaraLab/GWProt``

in a terminal window. 



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
   Examples/Example_2_Circular_Permutations
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