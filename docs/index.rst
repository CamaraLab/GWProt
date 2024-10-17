.. GWProt documentation master file


Welcome to GWProt's documentation!
===================================

**GWProt**  is a Python library for applying Gromov-Wasserstein to protein morphology developed by the `Cámara Lab <https://camara-lab.org/>`_



.. note::

   This project is under active development. Documentation and code will be updated continuously

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