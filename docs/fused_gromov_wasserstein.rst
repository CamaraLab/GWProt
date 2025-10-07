Fused Gromov-Wasserstein Correspondences
=========================================

Overview
--------
Fused Gromov-Wasserstein (FGW) [1]_ is a variant of the Gromov-Wasserstein distance that, 
in the context of GWProt, allows biochemical data to be incorporated into structural alignments.

In FGW, each residue is described not only by its spatial coordinates but also by additional 
data, such as hydrophobicity, isoelectric point, or a substitution score (e.g., from the BLOSUM62 matrix). 
FGW adds a penalty term for aligning biochemically dissimilar residues, enabling more 
biologically meaningful alignments.

Mathematical Formulation
------------------------
Let :math:`\delta(x, y)` denote the difference in the biochemical data associated with 
residues :math:`x` and :math:`y`. The *fused Gromov-Wasserstein distance* between 
proteins :math:`X` and :math:`Y` is defined as:

.. math::
   FGW(X, Y) = \min_T \frac{1}{2} \left( \sum_{i, j, k, l} \left[ \alpha |d_X(x_i, x_j) - d_Y(y_k, y_l)|^2 + (1 - \alpha) \delta(x_i, y_j) \right] T_{i, k} T_{j, l} \right)^{1/2}

where :math:`\alpha \in [0, 1]` determines the weight of the geometric cost relative 
to the biochemical penalty. Optimal values of :math:`\alpha` may vary depending on the 
type of biochemical data used.

References
----------
.. [1] Vayer, T., Chapel, L., Flamary, R., Tavenard, R., & Courty, N. (2019). Fused Gromov-Wasserstein distance for structured objects. Advances in Neural Information Processing Systems, 32.

