Distortion Scaling
==================

Overview
--------
When comparing the shapes of molecules, we are generally more interested in substructures that 
are geometrically close to one another, as these are more likely to interact chemically. This 
principle underlies the `template-modelling (TM) score <https://en.wikipedia.org/wiki/Template_modeling_score>`_ and 
`TM-align <https://zhanggroup.org/TM-align/TM-align>`_.

To implement this idea in GWProt, we modify the intra-protein distance matrices by applying a 
scaling function, causing the GW computation to give greater weight to nearby residues. This 
approach generally improves alignment accuracy without impacting runtime.

Scaling Function
----------------
We choose a *scaling function* :math:`f` such that :math:`f(0) = 0`, :math:`f` is strictly 
increasing, and :math:`f` is concave down. The square root function, :math:`f(x) = \sqrt{x}`, 
works well in practice. For each protein, we apply :math:`f` to all entries in its intra-protein 
distance matrix before running GW.

.. note::
   Distortion scaling can also be used with FGW in the same way as with GW. However, it is 
   recommended to use a larger value of ``alpha``.

References
----------
- Zhang, Y., & Skolnick, J. (2004). Scoring function for automated assessment of protein structure template quality. Proteins: Structure, Function, and Bioinformatics, 57(4), 702-710. (`TM-score paper <https://zhanggroup.org/TM-score/TMscore.pdf>`_)
- TM-align: https://zhanggroup.org/TM-align/




